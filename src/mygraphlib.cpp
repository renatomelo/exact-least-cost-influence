#include <math.h>
#include "myutils.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include<lemon/math.h>
#include "mygraphlib.h"
#include "myutils.h"
#include "mycolor.h"

#if __cplusplus >= 201103L
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif
using namespace std;

mygraphlibdefaultparameters mygraphlibparameters;

// The code below is divided in sections:
//     General Functions

void PulaBrancoComentario(ifstream &ifile)
{ char c;  string line;
  while (!ifile.eof()) {
    c = ifile.get();
    while ((c==' ') && (!ifile.eof())) c = ifile.get();
    if (!ifile.eof()) { if (c=='#')  getline(ifile,line);  else {ifile.unget(); break; } }
  }
}



bool WriteGraphGraphviz(Graph &g,
           NodeStringMap &vname, // vertex names
           EdgeStringMap &ename,  // edge names
           NodeColorMap &vcolor,   // vertex colors
           EdgeColorMap &acolor,    // edge colors
           string filename)
{
  ofstream out;
  string linha;

  out.open(filename.c_str());
  if (out.is_open()) return(false);

  out << "graph g {\n";
  out << "\tsize = \"8, 11\";\n";
  out << "\tnode [shape = \"circle\"];\n";
  for (NodeIt v(g); v!=INVALID; ++v) {
    linha = "\t";    linha += vname[v].c_str();   linha += " [color=";
    linha += ColorRGB(vcolor[v]);                linha += "];\n";
    out << linha;
  }
  for (EdgeIt a(g); a!=INVALID; ++a) {
    if (acolor[a]!=ColorCode("White")) {
      linha = "\t";
      linha += vname[g.u(a)].c_str() ;
      linha += "  -- ";
      linha += vname[g.v(a)].c_str();
      linha += " [label = \"";  linha += ename[a].c_str();  linha += "\"";
      linha += ", color=\""; linha += ColorRGB(acolor[a]); linha += "\" ];\n";
      out << linha;
    }
  }
  out << "}\n";
  out.close();
  return(true);
}




/* Dado um texto e um padrao, devolve a posicao no
   texto onde o padrao ocorre. Devolve (-1) se nao achar. */
int gr_busca_padrao(char *texto,char *busca)
{
  int i,j,tam_texto,tam_busca;
  char *p1,*p2;
  bool achou;
  tam_texto = strlen(texto);
  tam_busca = strlen(busca);
  for (i=0;i<tam_texto - tam_busca + 1;i++) {
      p1 = &texto[i];
      p2 = busca;
      for (j=0,achou = true;j<tam_busca && achou == true;j++)
    if (p1[j] != p2[j]) achou = false;
      if (achou) return(i);
  }
  return(-1);
}


Node GetNodeByName(Graph &g,
           NodeStringMap  &vname, // name of the nodes
           string vertexname)
{
  for (NodeIt v(g); v!=INVALID; ++v) {
    if (vname[v]==vertexname) return(v);
  }
  cout << "GetNodeByName: Error to obtain vertex node \"" << vertexname << "\"" << endl;
  exit(0);
}

// This routine uses a graphiz program to generate positions.
bool GenerateVertexPositions(Graph &g,
                 EdgeValueMap &custo,
                 NodePosMap   &posx,
                 NodePosMap   &posy)
{
  size_t t=0;
  double x,y;
  char tempname[1000],tempnamedot[1000],tempnameposdot[1000],cmd[1000];
  ofstream out;
  ifstream in;
  string linha,substring;

  (void) custo;// to avoid "non-used" parameter message.

  // obtain a temporary file name
  strcpy(tempname,".readgraphtempname");
  strcpy(tempnamedot,tempname);   strcat(tempnamedot,".dot");
  strcpy(tempnameposdot,tempname);strcat(tempnameposdot,"_pos.dot");

  out.open(tempnamedot);
  if (!out.is_open()) return(false);

  out << "graph g {\n";
  out << "\tsize = \"11, 11\";\n";
  out << "\tnode [shape = \"circle\"];\n";
  for (NodeIt v(g); v!=INVALID; ++v) {
    linha = "\t";    linha += IntToString(g.id(v));   linha += ";\n";
    out << linha;
  }
  for (EdgeIt a(g); a!=INVALID; ++a) {
    linha = "\t";  linha += IntToString(g.id(g.u(a)));
    linha += "  -- ";  linha += IntToString(g.id(g.v(a))); linha += ";\n";
    out << linha;
  }
  out << "}\n";
  out.close();
  sprintf(cmd,"%s -Goverlap=false %s -o %s",mygraphlibparameters.graphvizdrawingprogram.c_str(),tempnamedot,tempnameposdot);
  fflush(stdout);
  system(cmd); // gera outro arquivo do graphviz, mas com posições

  in.open(tempnameposdot);
  if (!in.is_open()) return(false);
  while (!in.eof()) {
    getline(in,linha);
    t = linha.find("{");
    if (t!=string::npos) break;
  }
  if (t==string::npos) {
    cout << "Temp Graphviz file is not appropriate for GenerateVertexPositions.\n";
    exit(0);
  }
  for (NodeIt v(g); v!=INVALID; ++v) {
    getline(in,linha);
    // avoid info on the graph, node or edges
    while ((!in.eof()) && ((linha.find("graph [")!=string::npos) ||
               (linha.find("node [")!=string::npos) ||
               (linha.find("edge [")!=string::npos) ||
               (linha.find(" -- ")!=string::npos))) {
      while ((!in.eof()) && (linha.find("];")==string::npos)) {
    string linha2;
    getline(in,linha2);
    linha += linha2;
      }
      getline(in,linha);
      while ((!in.eof()) && (linha.find("[")==string::npos)) {
    string linha2;
    getline(in,linha2);
    linha += linha2;
      }
    }
    if (linha.find("[")!=string::npos) {
      while (linha.find("];")==string::npos) {
    string linha2;
    getline(in,linha2);
    linha += linha2;
      }
    }
    t = linha.find("pos=\"");
    if (t!=string::npos) {
      stringstream s;
      int nodeid;
      sscanf(linha.c_str(),"%d",&nodeid);
      substring = linha.substr(t+5);
      sscanf(substring.c_str(),"%lf,%lf",&x,&y);
      //printf("%lf , %lf",x,y);
      for (NodeIt vv(g); vv!=INVALID; ++vv) {
    if (nodeid==g.id(vv)){
      posx[vv] = x;
      posy[vv] = y;
      //printf("interno: %d   ( %lf , %lf )\n",g.id(vv),posx[vv],posy[vv]);
      break;
    }
      }
    } else {
      printf("GenerateVertexPositions: Error to obtain vertex coordinates.\n");
      return(false);
    }
  }
  //for (NodeIt vv(g); vv!=INVALID; ++vv) {
  // printf(" %d   ( %lf , %lf )\n",g.id(vv),posx[vv],posy[vv]);
  //}
  //printf("-------------------------------------------------------\n");
  return(true);
}


// This routine uses a graphviz program to generate positions.
bool GenerateVertexPositions(Digraph &g,
                 ArcValueMap &custo,
                 DNodePosMap   &posx,
                 DNodePosMap   &posy)
{
  size_t t=0;
  double x,y;
  char tempname[1000],tempnamedot[1000],tempnameposdot[1000],cmd[1000];
  ofstream out;
  ifstream in;
  string linha,substring;

  (void) custo;// to avoid "non-used" parameter message.

  // obtain a temporary file name
  strcpy(tempname,".readgraphtempname");
  strcpy(tempnamedot,tempname);   strcat(tempnamedot,".dot");
  strcpy(tempnameposdot,tempname);strcat(tempnameposdot,"_pos.dot");

  out.open(tempnamedot);
  if (!out.is_open()) return(false);

  out << "digraph g {\n";
  out << "\tsize = \"11, 11\";\n";
  out << "\tnode [shape = \"circle\"];\n";
  for (DNodeIt v(g); v!=INVALID; ++v) {
    linha = "\t";    linha += IntToString(g.id(v));   linha += ";\n";
    out << linha;
  }
  for (ArcIt a(g); a!=INVALID; ++a) {
    linha = "\t";  linha += IntToString(g.id(g.source(a)));
    linha += "  -> ";  linha += IntToString(g.id(g.target(a))); linha += ";\n";
    out << linha;
  }
  out << "}\n";
  out.close();
  sprintf(cmd,"%s -Goverlap=false %s -o %s",mygraphlibparameters.graphvizdrawingprogram.c_str(),tempnamedot,tempnameposdot);
  fflush(stdout);
  system(cmd); // gera outro arquivo do graphviz, mas com posições

  in.open(tempnameposdot);
  if (!in.is_open()) return(false);
  while (!in.eof()) {
    getline(in,linha);
    t = linha.find("{");
    if (t!=string::npos) break;
  }
  if (t==string::npos) {
    cout << "Temp Graphviz file is not appropriate for GenerateVertexPositions.\n";
    exit(0);
  }
  for (DNodeIt v(g); v!=INVALID; ++v) {
    getline(in,linha);
    // avoid info on the graph, node or edges
    while ((!in.eof()) && ((linha.find("digraph [")!=string::npos) ||
               (linha.find("graph [")!=string::npos) ||
               (linha.find("node [")!=string::npos) ||
               (linha.find("edge [")!=string::npos) ||
               (linha.find(" -> ")!=string::npos))) {
      while ((!in.eof()) && (linha.find("];")==string::npos)) {
    string linha2;
    getline(in,linha2);
    linha += linha2;
      }
      getline(in,linha);
      while ((!in.eof()) && (linha.find("[")==string::npos)) {
    string linha2;
    getline(in,linha2);
    linha += linha2;
      }
    }
    if (linha.find("[")!=string::npos) {
      while (linha.find("];")==string::npos) {
    string linha2;
    getline(in,linha2);
    linha += linha2;
      }
    }
    t = linha.find("pos=\"");
    if (t!=string::npos) {
      stringstream s;
      int nodeid;
      sscanf(linha.c_str(),"%d",&nodeid);
      substring = linha.substr(t+5);
      sscanf(substring.c_str(),"%lf,%lf",&x,&y);
      // printf("teste %lf , %lf",x,y);
      for (DNodeIt vv(g); vv!=INVALID; ++vv) {
    if (nodeid==g.id(vv)){
      posx[vv] = x;
      posy[vv] = y;
      // printf("interno: %d   ( %lf , %lf )\n",g.id(vv),posx[vv],posy[vv]);
      break;
    }
      }
    } else {
      printf("GenerateVertexPositions: Error to obtain vertex coordinates.\n");
      return(false);
    }
  }
  //for (NodeIt vv(g); vv!=INVALID; ++vv) {
  // printf(" %d   ( %lf , %lf )\n",g.id(vv),posx[vv],posy[vv]);
  //}
  //printf("-------------------------------------------------------\n");
  return(true);
}


bool GenerateVertexPositions(Graph &g,NodePosMap  &px, NodePosMap &py){
  EdgeValueMap w(g);
  for (EdgeIt e(g); e!=INVALID; ++e) w[e] = 1.0;
  return(GenerateVertexPositions(g,w,px,py)); }

bool GenerateVertexPositions(Digraph &g,DNodePosMap  &px, DNodePosMap &py){
  ArcValueMap w(g);
  for (ArcIt a(g); a!=INVALID; ++a) w[a] = 1.0;
  return(GenerateVertexPositions(g,w,px,py)); }

bool GraphTable::GetNodeCoordinates(string colx, NodePosMap &posx, string coly, NodePosMap &posy){
    // read from file
    if((this->NodeTable->columnindex(colx) != -1) && (this->NodeTable->columnindex(coly) != -1)){
        this->GetColumn(colx,posx);
        this->GetColumn(coly,posy);
    }
    // no positions in the file, call a routine to create positions for the nodes
    else
        GenerateVertexPositions(this->g, posx, posy);

    return true;
}

bool GetNodeCoordinates(DigraphTable &GT, string colx, DNodePosMap &posx, string coly, DNodePosMap &posy)
{
  if ((GT.NodeTable->columnindex(colx)!=-1)&&
      (GT.NodeTable->columnindex(coly)!=-1)){
    GT.GetColumn(colx,posx);
    GT.GetColumn(coly,posy);
  } else GenerateVertexPositions(GT.g,posx,posy);
  return(true);
}

// Trocar depois para Templates para serem uma rotina, jah que soh muda o tipo

bool InitGraphTable(StringTable &T,Graph &g,LineToNodeMap &nodemap)
{  // Make a correspondence between row i and node v
  nodemap.reserve(T.nrows); nodemap.resize(T.nrows);
  for (int i=0;i<T.nrows;i++){ Node v=g.addNode();  nodemap[i] = v; } return(true);}

bool InitGraphTable(StringTable &T,Graph &g,StringToNodeMap &smap,
            string colname_v1,string colname_v2,LineToEdgeMap &edgemap)
{  // Make a correspondence between row i and edge e
  edgemap.reserve(T.nrows); edgemap.resize(T.nrows);
  int col_u=T.columnindex(colname_v1),col_v=T.columnindex(colname_v2);
  string name_u,name_v;
  for (int i=0;i<T.nrows;i++){
    T.entry(i,col_u,name_u);  T.entry(i,col_v,name_v);
    Node u=smap[name_u],v=smap[name_v];
    Edge e=g.addEdge(u,v); edgemap[i] = e; } return(true);}

// Insert nodes in the graph
bool InitGraphTable(StringTable &T,Digraph &g,LineToDNodeMap &dnodemap)
{  // Make a correspondence between row i and node v
  dnodemap.reserve(T.nrows); dnodemap.resize(T.nrows);
  for (int i=0;i<T.nrows;i++){ DNode v=g.addNode();  dnodemap[i] = v; } return(true);}

// insert arcs in the graph
bool InitGraphTable(StringTable &T,Digraph &g,StringToDNodeMap &smap,
            string colname_source,string colname_target,LineToArcMap &arcmap)
{  // Make a correspondence between row i and arc a
  arcmap.reserve(T.nrows); arcmap.resize(T.nrows);
  int col_source=T.columnindex(colname_source),col_target=T.columnindex(colname_target);
  //cout << "nrows = " << T.nrows << endl;
  //cout << "source and target: " << col_source << "   " << col_target << endl;
  string name_source,name_target;
  for (int i=0;i<T.nrows;i++){
    T.entry(i,col_source,name_source);  T.entry(i,col_target,name_target);
    //cout << "namesource and nametarget: " << name_source << "   " << name_target << endl;
    DNode source=smap[name_source],target=smap[name_target];
    //cout << "sourceid and targetid: " << g.id(source) << "   " << g.id(target) << endl;
    Arc a=g.addArc(source,target); arcmap[i] = a; } return(true);}

//------------------------------------------------------------
// Read a column with Node info
bool ReadGraphColumn(StringTable &T,LineToNodeMap &nodemap,string colname,NodeValueMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[nodemap[i]]); return(true);}

bool ReadGraphColumn(StringTable &T,LineToNodeMap &nodemap,string colname,NodeStringMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[nodemap[i]]); return(true);}

bool ReadGraphColumn(StringTable &T,LineToNodeMap &nodemap,string colname,NodeIntMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[nodemap[i]]); return(true);}
//------------------------------------------------------------
// Read a column with Edge info
bool ReadGraphColumn(StringTable &T,LineToEdgeMap &edgemap,string colname,EdgeValueMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[edgemap[i]]); return(true);}
bool ReadGraphColumn(StringTable &T,LineToEdgeMap &edgemap,string colname,EdgeStringMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[edgemap[i]]); return(true);}
bool ReadGraphColumn(StringTable &T,LineToEdgeMap &edgemap,string colname,EdgeIntMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[edgemap[i]]); return(true);}
//------------------------------------------------------------
// Read a column with DNode info
bool ReadGraphColumn(StringTable &T,LineToDNodeMap &nodemap,string colname,DNodeValueMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[nodemap[i]]); return(true);}
bool ReadGraphColumn(StringTable &T,LineToDNodeMap &nodemap,string colname,DNodeStringMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[nodemap[i]]); return(true);}
bool ReadGraphColumn(StringTable &T,LineToDNodeMap &nodemap,string colname,DNodeIntMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[nodemap[i]]); return(true);}
//------------------------------------------------------------
// Read a column with Arc info
bool ReadGraphColumn(StringTable &T,LineToArcMap &arcmap,string colname,ArcValueMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[arcmap[i]]); return(true);}
bool ReadGraphColumn(StringTable &T,LineToArcMap &arcmap,string colname,ArcStringMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[arcmap[i]]); return(true);}
bool ReadGraphColumn(StringTable &T,LineToArcMap &arcmap,string colname,ArcIntMap &read) {
  int col = T.columnindex(colname);
  for (int i=0;i<T.nrows;i++) T.entry(i,col,read[arcmap[i]]); return(true);}
//------------------------------------------------------------

bool ReadGraph(string filename,
              Graph &g,
              NodeStringMap& vname,
              NodePosMap& posx,
              NodePosMap& posy,
              EdgeValueMap& weight)
{
  string type = GetGraphFileType(filename);
  if (type!="graph"){cout<<"Error: Unknown type of graph: "<<type<<endl;exit(1);}
  GraphTable GT(filename,g); // Read the graph (only nodes and edges)
  bool ok = GT.GetColumn(mygraphlibparameters.inputgraphtable_nodename,vname);
  ok = ok && GT.GetColumn("weight",weight);
  ok = ok && GT.GetNodeCoordinates(mygraphlibparameters.inputgraphtable_posx,posx,
                   mygraphlibparameters.inputgraphtable_posy,posy);
  return(ok);
}

bool ReadGraph(string filename,
           Graph &g,
           NodeStringMap& vname,
           NodePosMap& posx,
           NodePosMap& posy,
           NodeValueMap& weight)
{
  string type = GetGraphFileType(filename);
  if (type!="graph"){cout<<"Error: Unknown type of graph: "<<type<<endl;exit(1);}
  GraphTable GT(filename,g); // Read the graph (only nodes and edges)
  bool ok = GT.GetColumn(mygraphlibparameters.inputgraphtable_nodename,vname);
  ok = ok && GT.GetColumn("weight",weight);
  ok = ok && GT.GetNodeCoordinates(mygraphlibparameters.inputgraphtable_posx,posx,
                   mygraphlibparameters.inputgraphtable_posy,posy);
  return(ok);
}


bool ReadDigraph(string filename,
         Digraph &g,
         DNodeStringMap& vname,
         DNodePosMap& posx,
         DNodePosMap& posy,
         ArcValueMap& weight)
{
  string type = GetGraphFileType(filename);
  if (type!="digraph"){cout<<"Error: Unknown type of digraph: "<<type<<endl;exit(1);}
  DigraphTable GT(filename,g); // Read the graph (only nodes and edges)
  bool ok = GT.GetColumn("nodename",vname);
  ok = ok && GT.GetColumn("weight",weight);
  ok = ok && GetNodeCoordinates(GT,"posx",posx,"posy",posy);
  return(ok);
}



//Generate a random complete euclidean Graph
bool GenerateRandomEuclideanDigraph(Digraph &g,
              DNodeStringMap &vname, // node name
              DNodePosMap    &px, // x-position of the node
              DNodePosMap    &py, // y-position of the node
              ArcValueMap  & weight, // weight of edges
              int n, // number of nodes
              double SizeX, // coordinate x is a random number in [0,SizeX)
              double SizeY) // coordinate y is a random number in [0,SizeY)
{
  DNode *V;
  V = new DNode[n];
  if (V==NULL){
    cout << "Memory allocation error, number of nodes " << n << " too large\n";
    exit(0);}

  for (int i=0;i<n;i++) {   // insert nodes (random points in [0,100] x [0,100] )
    V[i] = g.addNode();    // new node
    px[V[i]] = SizeX*drand48();
    py[V[i]] = SizeY*drand48();
    vname[V[i]] = IntToString(i+1);// name of node is i+1
  }
  for (int i=0;i<n;i++)
    for (int j=i+1;j<n;j++) {
      Arc e = g.addArc(V[i],V[j]);  // generate new arc from v_i to v_j
      weight[e] = sqrt(pow(px[V[i]]-px[V[j]],2) + pow(py[V[i]]-py[V[j]],2));
      Arc f = g.addArc(V[j],V[i]);  // generate new arc from v_j to v_i
      weight[f] = weight[e];
    }
  delete[] V;
  return(true);
}


// Generate a triangulated Graph, building the Delaunay
// triangulation of random points
// Uses the geompack program, available in
// http://people.sc.fsu.edu/~jburkardt/cpp_src/geompack/geompack.html
bool GenerateTriangulatedGraph(Graph &g, // return with generated graph
              NodeStringMap &vname, // return with name of the nodes
              NodePosMap& px, // return with x-position of the nodes
              NodePosMap& py, // return with y-position of the nodes
              EdgeValueMap& weight, // return with weight of edges
              int n, // number of nodes
              double SizeX, // coordinate x is a random number in [0,SizeX)
              double SizeY) // coordinate y is a random number in [0,SizeY)
{
  int i; // n=number of nodes
  int ntri; // number of Delaunay triangles
  Node *V = new Node[n];
  double *p = new double[2*n+2];// node coodinates are (x;y) = ( p[2*i] ; p[2*i+1] )
  int *tri = new int[6*n]; // Each 3 elements are the indexes of a triangle
  int *tri_nabe = new int[6*n];

  if ((V==NULL)||(p==NULL)||(tri==NULL)||(tri_nabe==NULL)){
    cout << "Memory allocation error, number of nodes " << n << " too large\n";
    exit(0);}

  for (i=0;i<n;i++) {
    V[i] = g.addNode();    // gera um vértice nó do grafo
    px[V[i]] = SizeX*drand48();  // nodes are random points
    py[V[i]] = SizeY*drand48();
    p[2*i]   = px[V[i]];
    p[2*i+1] = py[V[i]];
    vname[V[i]] = IntToString(i+1);  // name of the node is i+1
  }
  if (r8tris2 ( n, p, &ntri, tri,  tri_nabe )) {  printf("ERROR\n");Pause(); }
  for (i=0;i<ntri;i++) {
    int a,b,c;
    a = tri[3*i]-1; b = tri[3*i+1]-1; c = tri[3*i+2]-1;
    // each triangle if formed with nodes  V[a] , V[b] , V[c]
    // insert edges without duplications
    if ((findEdge(g,V[a],V[b])==INVALID) && (findEdge(g,V[b],V[a])==INVALID)){
      Edge e = g.addEdge(V[a],V[b]);
      weight[e] = sqrt(pow(px[V[a]]-px[V[b]],2) + pow(py[V[a]]-py[V[b]],2));
    }
    if ((findEdge(g,V[a],V[c])==INVALID)  && (findEdge(g,V[c],V[a])==INVALID)){
      Edge e = g.addEdge(V[a],V[c]);
      weight[e] = sqrt(pow(px[V[a]]-px[V[c]],2) + pow(py[V[a]]-py[V[c]],2));
    }
    if ((findEdge(g,V[b],V[c])==INVALID) && (findEdge(g,V[c],V[b])==INVALID)) {
      Edge e = g.addEdge(V[b],V[c]);
      weight[e] = sqrt(pow(px[V[b]]-px[V[c]],2) + pow(py[V[b]]-py[V[c]],2));
    }
  }
  delete[] V;
  delete[] p;
  delete[] tri;
  delete[] tri_nabe;
  return(true);
}

// Generate a triangulated Digraph, building the Delaunay
// triangulation of random points. Each edge of the Delaunay triangulation
// leads to two arcs (in both senses)
// Uses the geompack program, available in
// http://people.sc.fsu.edu/~jburkardt/cpp_src/geompack/geompack.html
bool GenerateTriangulatedDigraph(Digraph &g,
         DNodeStringMap &vname, // name of the nodes
         DNodePosMap  &px, // x-position of the nodes
         DNodePosMap  &py, // y-position of the nodes
         ArcValueMap   & weight, // weight of edges
         int n, // number of nodes
         double SizeX, // coordinate x is a random number in [0,SizeX)
         double SizeY) // coordinate y is a random number in [0,SizeY)
{
 int i; // n=number of nodes
 DNode *V = new DNode[n];
 double *p = new double[2*n+2];// node coodinates are (x;y) = ( p[2*i] ; p[2*i+1] )
 int *tri = new int[6*n]; // Each 3 elements are the indexes of a triangle
 int ntri;
 int *tri_nabe = new int[6*n];

  if ((V==NULL)||(p==NULL)||(tri==NULL)||(tri_nabe==NULL)){
    cout << "Memory allocation error, number of nodes " << n << " too large\n";
    exit(0);}

 for (i=0;i<n;i++) {   // insere os vértices (pontos aleatórios no plano [0,SizeX]x[0,SizeY] )
    V[i] = g.addNode(); // gera um vértice nó do grafo
    px[V[i]] = SizeX*drand48();  // nodes are random points
    py[V[i]] = SizeY*drand48();
    p[2*i]   = px[V[i]];
    p[2*i+1] = py[V[i]];
    vname[V[i]] = IntToString(i+1);// name of node is i+1
  }
  if (r8tris2 ( n, p, &ntri, tri,  tri_nabe )) {  printf("ERROR\n");Pause(); }
  for (i=0;i<ntri;i++) {
    int a,b,c;
    a = tri[3*i]-1; b = tri[3*i+1]-1; c = tri[3*i+2]-1;
    // each triangle if formed with nodes  V[a] , V[b] , V[c]
    // insert arcs without duplications in both senses
    if (findArc(g,V[a],V[b])==INVALID){
      Arc e = g.addArc(V[a],V[b]);
      weight[e] = sqrt(pow(px[V[a]]-px[V[b]],2) + pow(py[V[a]]-py[V[b]],2));
    }
    if (findArc(g,V[b],V[a])==INVALID){
      Arc e = g.addArc(V[b],V[a]);
      weight[e] = sqrt(pow(px[V[b]]-px[V[a]],2) + pow(py[V[b]]-py[V[a]],2));
    }
    if (findArc(g,V[a],V[c])==INVALID){
      Arc e = g.addArc(V[a],V[c]);
      weight[e] = sqrt(pow(px[V[a]]-px[V[c]],2) + pow(py[V[a]]-py[V[c]],2));
    }
    if (findArc(g,V[c],V[a])==INVALID){
      Arc e = g.addArc(V[c],V[a]);
      weight[e] = sqrt(pow(px[V[c]]-px[V[a]],2) + pow(py[V[c]]-py[V[a]],2));
    }
    if (findArc(g,V[b],V[c])==INVALID) {
      Arc e = g.addArc(V[b],V[c]);
      weight[e] = sqrt(pow(px[V[b]]-px[V[c]],2) + pow(py[V[b]]-py[V[c]],2));
    }
    if (findArc(g,V[c],V[b])==INVALID) {
      Arc e = g.addArc(V[c],V[b]);
      weight[e] = sqrt(pow(px[V[c]]-px[V[b]],2) + pow(py[V[c]]-py[V[b]],2));
    }
  }
  delete[] V;
  delete[] p;
  delete[] tri;
  delete[] tri_nabe;
  return(true);
}


//Generate a random complete euclidean Graph
bool GenerateRandomEuclideanGraph(Graph &g,
          NodeStringMap &vname, // node names
          NodePosMap& px, // x-position of the nodes
          NodePosMap& py, // y-position of the nodes
          EdgeValueMap& weight, // weight of edges
          int n, // number of nodes
          double SizeX, // coordinate x is a random number in [0,SizeX)
          double SizeY) // coordinate y is a random number in [0,SizeY)
{
  int i,j; // n=number of nodes
  Node *V;
  V = new Node[n];
  if (V==NULL){
    cout << "Memory allocation error, number of nodes " << n << " too large\n";
    exit(0);}

  for (i=0;i<n;i++) {   // insert nodes (random points in [0,100] x [0,100] )
    V[i] = g.addNode(); // generate a new node
    px[V[i]] = SizeX*drand48();
    py[V[i]] = SizeY*drand48();
    vname[V[i]] = IntToString(i+1); // name of the node is i+1
  }
  for (i=0;i<n;i++)
    for (j=i+1;j<n;j++) {
      Edge e = g.addEdge(V[i],V[j]);  // generate an edge
      weight[e] = sqrt(pow(px[V[i]]-px[V[j]],2) + pow(py[V[i]]-py[V[j]],2));
    }
  delete[] V;
  return(true);
}







// obtain a mincut separating vertices 's' and 't' for undirected graphs
// The input is given by the graph 'g'. Each edge of g has a "weight".
// The returned cut is given by the vector of nodes 'cut' (boolean
// vector, nodes in one side have value false and nodes in the other side
// have value true).
double MinCut(Graph &g, EdgeValueMap &weight, Node &s,Node &t, CutMap &cut)
{
  Preflow<Graph, EdgeValueMap> pf(g, weight, s, t);
  pf.runMinCut();
  pf.minCutMap(cut);
  return (pf.flowValue());
}

// Obtain a mininum cut for directed graphs from s to t.
// The returned cut is given by the vector of nodes 'cut' (boolean
// vector: nodes v in the same side of s have cut[v]=true, otherwise cut[v]=false.
double DiMinCut(Digraph &g, ArcValueMap &weight, DNode &s,DNode &t, DCutMap &vcut)
{
  Preflow<Digraph, ArcValueMap> preflow_test(g, weight, s, t);
  preflow_test.run();
  preflow_test.minCutMap(vcut);
  return (preflow_test.flowValue());
}




void PrintGraph(Graph &g, NodeStringMap &vname, EdgeValueMap &graphweight)
{
  int Nnodes = countNodes(g); // number of nodes in the input graph
  int Nedges = countEdges(g); // number of edges in the input graph

  printf("-------------------------------------------------------\n");
  printf("Number of nodes: %d\n",Nnodes);
  printf("Number of edges: %d\n",Nedges);
  for (NodeIt v(g); v!=INVALID; ++v) printf("%s\n",vname[v].c_str());
  printf("\n");
  printf("-------------------------------------------------------\n");
  for (EdgeIt a(g); a!=INVALID; ++a)
    printf("%s -- %s  %lf\n",vname[g.u(a)].c_str(),vname[g.v(a)].c_str(), graphweight[a]);
  printf("\n");
}


GraphAttributes::GraphAttributes(Graph &graph, NodeStringMap &vname,
                 NodePosMap &posx, NodePosMap &posy):
  g(graph),node_name(graph),px(graph),py(graph),node_attrib(graph),edge_attrib(graph)
{
  factor = 1.0;
  graph_attrib = "";
  for (NodeIt v(this->g);v!=INVALID;++v) {
    node_name[v] = vname[v];
    px[v] = posx[v];
    py[v] = posy[v];
    node_attrib[v] = "";  }
  for (EdgeIt e(this->g); e != INVALID; ++e) edge_attrib[e] = "";
}

GraphAttributes::GraphAttributes(GraphTable &GT):
  g(GT.g),node_name(GT.g),px(GT.g),py(GT.g),node_attrib(GT.g),edge_attrib(GT.g)
{
  factor = 1.0;
  graph_attrib = "";
  GT.GetColumn(mygraphlibparameters.inputgraphtable_nodename,node_name);
  GT.GetNodeCoordinates(mygraphlibparameters.inputgraphtable_posx,px,
                mygraphlibparameters.inputgraphtable_posy,py);
  for (NodeIt v(this->g);v!=INVALID;++v) node_attrib[v] = "";
  for (EdgeIt e(this->g); e != INVALID; ++e) edge_attrib[e] = "";

  GT.GetColumn("nodename",this->node_name);
  GT.GetNodeCoordinates(mygraphlibparameters.inputgraphtable_posx,this->px,
                mygraphlibparameters.inputgraphtable_posy,this->py);
}


// This routine visualize a graph using a pdf viewer. It uses a graphviz (from
// graphviz.org) to generate a pdf file and a program to view the pdf file. The
// pdf viewer name is given in the viewername parameter.
bool GraphAttributes::View() // text displayed below the figure
{
  char tempname[1000],cmd[1000],outputname[1000];
  FILE *fp;
  string str;
  double minpx=DBL_MAX,minpy=DBL_MAX,maxpx=-DBL_MAX,maxpy=-DBL_MAX,delta,normalfactor;

  // obtain a temporary file name
  strcpy(tempname,".viewgraphtempname");sprintf(outputname,"%s.pdf",tempname);
  fp = fopen(tempname,"w+");
  if (fp==NULL) {cout << "Error to open temporary file to visualize graph.\n"; return(false);}
  for (NodeIt v(this->g); v!=INVALID; ++v) {
    if (this->px[v] < minpx) minpx = this->px[v];
    if (this->px[v] > maxpx) maxpx = this->px[v];
    if (this->py[v] < minpy) minpy = this->py[v];
    if (this->py[v] > maxpy) maxpy = this->py[v];
  }
  normalfactor = 20/this->factor; // quanto maior, menor o desenho do vértice
  delta = fmax(maxpx - minpx,maxpy - minpy);
  // Generate a text file with the graph format of graphviz program
  fprintf(fp,"graph g {\n");
  fprintf(fp,"\t %s \n", this->graph_attrib.c_str()); // put graph attributes
  fprintf(fp,"\tnode [ %s ];\n", this->default_node_attrib.c_str()); // put default node attributes
  fprintf(fp,"\tedge [ %s ];\n", this->default_edge_attrib.c_str()); // put default edge attributes
  for (NodeIt v(this->g); v!=INVALID; ++v) {
    fprintf(fp,"\t%s [ ",this->node_name[v].c_str());
    fprintf(fp,"%s ",this->node_attrib[v].c_str());
    fprintf(fp,"pos=\"%lf,%lf!\" ",normalfactor*(this->px[v]-minpx)/delta,normalfactor*(this->py[v]-minpy)/delta);
    fprintf(fp," ];\n");
  }
  for (EdgeIt e(this->g); e!=INVALID; ++e) {
    if (IsSuffix(this->edge_attrib[e],"color=\"Invis\"")) continue;// only to speed graphviz exec.
    fprintf(fp,"\t%s  -- %s [",this->node_name[this->g.u(e)].c_str(),this->node_name[this->g.v(e)].c_str());
    fprintf(fp,"%s ",this->edge_attrib[e].c_str());
    fprintf(fp," ];\n");
  }

  fprintf(fp,"}\n");
  fclose(fp);
  sprintf(cmd,"%s -Tpdf %s -o %s",mygraphlibparameters.graphvizdrawingprogram.c_str(),tempname,outputname);
system(cmd);
  //cout << "Grafo em "<< tempname << "\n";
  view_pdffile(outputname);
  //pause();
  return(true);
}

// This routine visualize a graph using a pdf viewer. It uses a graphviz (from
// graphviz.org) to generate a pdf file and a program to view the pdf file. The
// pdf viewer name is given in the viewername parameter.
bool GraphAttributes::Save(string name) // text displayed below the figure
{
  char tempname[1000],cmd[1000],outputname[1000];
  FILE *fp;
  string str;
  double minpx=DBL_MAX,minpy=DBL_MAX,maxpx=-DBL_MAX,maxpy=-DBL_MAX,delta,normalfactor;

  // obtain a temporary file name
  strcpy(tempname,name.c_str());sprintf(outputname,"%s.pdf",tempname);
  fp = fopen(tempname,"w+");
  if (fp==NULL) {cout << "Error to open temporary file to visualize graph.\n"; return(false);}
  for (NodeIt v(this->g); v!=INVALID; ++v) {
    if (this->px[v] < minpx) minpx = this->px[v];
    if (this->px[v] > maxpx) maxpx = this->px[v];
    if (this->py[v] < minpy) minpy = this->py[v];
    if (this->py[v] > maxpy) maxpy = this->py[v];
  }
  normalfactor = 20/this->factor; // quanto maior, menor o desenho do vértice
  delta = fmax(maxpx - minpx,maxpy - minpy);
  // Generate a text file with the graph format of graphviz program
  fprintf(fp,"graph g {\n");
  fprintf(fp,"\t %s \n", this->graph_attrib.c_str()); // put graph attributes
  fprintf(fp,"\tnode [ %s ];\n", this->default_node_attrib.c_str()); // put default node attributes
  fprintf(fp,"\tedge [ %s ];\n", this->default_edge_attrib.c_str()); // put default edge attributes
  for (NodeIt v(this->g); v!=INVALID; ++v) {
    fprintf(fp,"\t%s [ ",this->node_name[v].c_str());
    fprintf(fp,"%s ",this->node_attrib[v].c_str());
    fprintf(fp,"pos=\"%lf,%lf!\" ",normalfactor*(this->px[v]-minpx)/delta,normalfactor*(this->py[v]-minpy)/delta);
    fprintf(fp," ];\n");
  }
  for (EdgeIt e(this->g); e!=INVALID; ++e) {
    if (IsSuffix(this->edge_attrib[e],"color=\"Invis\"")) continue;// only to speed graphviz exec.
    fprintf(fp,"\t%s  -- %s [",this->node_name[this->g.u(e)].c_str(),this->node_name[this->g.v(e)].c_str());
    fprintf(fp,"%s ",this->edge_attrib[e].c_str());
    fprintf(fp," ];\n");
  }

  fprintf(fp,"}\n");
  fclose(fp);
  sprintf(cmd,"%s -Tpdf %s -o %s",mygraphlibparameters.graphvizdrawingprogram.c_str(),tempname,outputname);
  system(cmd);
  return(true);
}

DigraphAttributes::DigraphAttributes(Digraph &digraph, DNodeStringMap &vname,
                 DNodePosMap &posx, DNodePosMap &posy):
g(digraph),node_name(vname),px(posx),py(posy),node_attrib(digraph),arc_attrib(digraph)
{
  factor = 1.0;
  digraph_attrib = "";
  for (DNodeIt v(this->g);v!=INVALID;++v) node_attrib[v] = "";
  for (ArcIt e(this->g); e != INVALID; ++e) arc_attrib[e] = "";
}


// This routine visualize a graph using a pdf viewer. It uses a graphviz (from
// graphviz.org) to generate a pdf file and a program to view the pdf file. The
// pdf viewer name is given in the viewername parameter.
bool DigraphAttributes::View() // text displayed below the figure
{
  char tempname[1000],cmd[1000],outputname[1000];
  FILE *fp;
  string str;
  double minpx=DBL_MAX,minpy=DBL_MAX,maxpx=-DBL_MAX,maxpy=-DBL_MAX,delta,normalfactor;

  // obtain a temporary file name
  strcpy(tempname,".viewdigraphtempname");
  sprintf(outputname,"%s.pdf",tempname);
  fp = fopen(tempname,"w+");
  if (fp==NULL) {cout << "Error to open temporary file to visualize digraph.\n"; return(false);}
  for (DNodeIt v(this->g); v!=INVALID; ++v) {
    if (this->px[v] < minpx) minpx = this->px[v];
    if (this->px[v] > maxpx) maxpx = this->px[v];
    if (this->py[v] < minpy) minpy = this->py[v];
    if (this->py[v] > maxpy) maxpy = this->py[v];
  }
  normalfactor = 20/this->factor; // quanto maior, menor o desenho do vértice
  delta = fmax(maxpx - minpx,maxpy - minpy);
  // Generate a text file with the graph format of graphviz program
  fprintf(fp,"digraph g {\n");
  //fprintf(fp,"\tsize = \"10, 10\";\n");
  fprintf(fp,"\t %s \n", this->digraph_attrib.c_str()); // put graph attributes
  fprintf(fp,"\tnode [ %s ];\n", this->default_node_attrib.c_str()); // put default node attributes
  // For digraphs, the graphviz use also the term "edge" in the file format
  fprintf(fp,"\tedge [ %s ];\n", this->default_arc_attrib.c_str()); // put default arc attributes
  for (DNodeIt v(this->g); v!=INVALID; ++v) {
    fprintf(fp,"\t%s [ ",this->node_name[v].c_str());
    fprintf(fp,"%s ",this->node_attrib[v].c_str());
    fprintf(fp,"pos=\"%lf,%lf!\" ",normalfactor*(this->px[v]-minpx)/delta,normalfactor*(this->py[v]-minpy)/delta);
    fprintf(fp," ];\n");
  }
  for (ArcIt e(this->g); e!=INVALID; ++e) {
    fprintf(fp,"\t%s  -> %s [",this->node_name[this->g.source(e)].c_str(),this->node_name[this->g.target(e)].c_str());
    fprintf(fp,"%s ",this->arc_attrib[e].c_str());
    fprintf(fp," ];\n");
  }

  fprintf(fp,"}\n");
  fclose(fp);
  sprintf(cmd,"%s -Tpdf %s -o %s",mygraphlibparameters.graphvizdrawingprogram.c_str(),tempname,outputname);
system(cmd);
  //cout << "Grafo em "<< tempname << "\n";
  view_pdffile(outputname);
  //pause();
  return(true);
}


bool randomgraph(Graph &g,NodeStringMap &vname, int n, double p)
{
  vector<Node> v(n);
  for (int i=0;i<n;i++){ v[i] = g.addNode(); vname[v[i]] = IntToString(i+1);  }
  for (int i=0;i<n-1;i++) for (int j=i+1;j<n;j++) if (drand48()<p) g.addEdge(v[i],v[j]);
  return(true);
}

bool randomdigraph(Digraph &g,DNodeStringMap &vname, int n, double p)
{
  vector<DNode> v(n);
  for (int i=0;i<n;i++){ v[i] = g.addNode(); vname[v[i]] = IntToString(i+1);  }
  for (int i=0;i<n-1;i++)
    for (int j=i+1;j<n;j++)
      if (drand48()<p) {
    // each direction with half of the probability p
    if (drand48() < 1.0/2) g.addArc(v[i],v[j]); else g.addArc(v[j],v[i]); }

  return(true);
}

// ================================================================
int globalcounter=0;

void ADJMAT_FreeNotNull(void *p){  if (p) free(p);  }

// Define an adjacency matrix, so as we have a fast access for the edges of a graph,
// given a pair of vertices. This is mainly used in the subroutine 2opt, for
// dense graphs.
// The adjacency matrix is stored in a strictly triangular inferior matrix.

AdjacencyMatrix::AdjacencyMatrix(Graph &graph,EdgeValueMap &graphweight,double NonEdgValue):
  Node2Index(graph),Edge2Index(graph)
{
  int i;
  g = &graph;
  NonEdgeValue = NonEdgValue;
  weight = &graphweight;
  Nnodes = countNodes(graph); // number of nodes in the input graph
  Nedges = countEdges(graph); // number of edges in the input graph
  Nmatrix = (Nnodes*(Nnodes-1))/2; // no. of edges/elem. in strict. triang. inf. matrix

  AdjMatrix = (double *) malloc(sizeof(double)*Nmatrix);
  Index2Node = (Node *) malloc(sizeof(Node)*Nnodes);
  Index2Edge = (Edge *) malloc(sizeof(Edge)*Nedges);

  if ((AdjMatrix==NULL)||(Index2Node==NULL)||(Index2Edge==NULL)) {
    cout << "Out of memory in constructor of AdjacencyMatrix\n";
    ADJMAT_FreeNotNull(AdjMatrix); ADJMAT_FreeNotNull(Index2Node); ADJMAT_FreeNotNull(Index2Edge);
    exit(0);}

  i = 0;
  for (NodeIt v(*g); v != INVALID; ++v) {
    Index2Node[i] = v;
    AdjacencyMatrix::Node2Index[v]=i;
    i++;
  }

  // Initially all edges have infinity weight
  for (int i=0;i<Nmatrix;i++) AdjMatrix[i] = NonEdgeValue;
  // Then, update the existing edges with the correct weight
  for (EdgeIt e(graph); e != INVALID; ++e) {
    Node u,v;    int i_u,i_v;
    u = graph.u(e);  v = graph.v(e);  // obtain the extremities of e
    i_u = Node2Index[u];
    i_v = Node2Index[v];
    if (i_u > i_v) AdjMatrix[i_u*(i_u-1)/2+i_v] = graphweight[e];
    else if (i_u < i_v) AdjMatrix[i_v*(i_v-1)/2+i_u] = graphweight[e];
    else {
      cout << "Out of memory in constructor of AdjacencyMatrix\n";
      exit(0);}
  }
}

double AdjacencyMatrix::Cost(Node u,Node v)
{
  int i_u,i_v;
  i_u = Node2Index[u];
  i_v = Node2Index[v];
  globalcounter++;
  //cout << "globalcounter1 = " << globalcounter << endl;
  try{
    if (i_u > i_v) return(AdjMatrix[i_u*(i_u-1)/2+i_v]);
    else if (i_u < i_v) return(AdjMatrix[i_v*(i_v-1)/2+i_u]);
    else return(NonEdgeValue);
  }catch (...) {
    cout << "Exception: " << globalcounter << endl;
    exit(1);
  }
}

double AdjacencyMatrix::Cost(Edge e)
{
  int i_u,i_v;
  Node u,v;
  globalcounter++;
  //cout << "globalcounter2 = " << globalcounter << endl;
  u = (*g).u(e);  v = (*g).v(e);
  i_u = Node2Index[u];
  i_v = Node2Index[v];
  if (i_u > i_v) return(AdjMatrix[i_u*(i_u-1)/2+i_v]);
  else return(AdjMatrix[i_v*(i_v-1)/2+i_u]);
}


AdjacencyMatrix::~AdjacencyMatrix()
{
  ADJMAT_FreeNotNull(AdjMatrix); ADJMAT_FreeNotNull(Index2Node); ADJMAT_FreeNotNull(Index2Edge);
}



// Generate a digraph from a given graph.
// For each edge e={u,v} in the graph, the digraph will have two arcs,
// (u,v), given in edgemap1[e], and (v,u) given by edgemap2[e]
bool Graph2Digraph(Graph &g, Digraph &dg, NodeDNodeMap &nodemap, EdgeArcMap &edgemap1, EdgeArcMap &edgemap2) {
  for (NodeIt v(g);v!=INVALID;++v) nodemap[v] = dg.addNode();
  for (EdgeIt e(g);e!=INVALID;++e) {
    edgemap1[e] = dg.addArc(nodemap[g.u(e)],nodemap[g.v(e)]);
    edgemap2[e] = dg.addArc(nodemap[g.v(e)],nodemap[g.u(e)]); }
  return(true);
}

// Heuristic to proper color the nodes of a graph.
// It is guaranteed that any planar graph is colored by at most 6 colors.
// ncolors returns with the number of used colors.
// color[v] returns with a color in {0,...,ncolors-1}
bool GreedyColoring(Graph &g,int &ncolors,NodeIntMap &nodecolor)
{
  int n=countNodes(g),maxdg;
  NodeNodeMap previous(g);
  NodeNodeMap next(g);
  NodeIntMap degree(g);
  NodeBoolMap removed(g);
  for (NodeIt v(g);v!=INVALID;++v){
    removed[v]=false; degree[v]=0; previous[v]=INVALID;  next[v]=INVALID; }
  for (EdgeIt e(g);e!=INVALID;++e){degree[g.u(e)]++; degree[g.v(e)]++;}
  maxdg=Max(g,degree);
  maxdg++; // all degrees are in {0,...,maxdg-1}
  vector<Node>listdegree(maxdg); // listdegree[d]: circular list of nodes with degree d
  for (int i=0;i<maxdg;i++) listdegree[i]=INVALID;
  for (NodeIt v(g);v!=INVALID;++v) {
    int d = degree[v]; // insert v in the list listdegree[d]
    if (listdegree[d]==INVALID){ next[v]=v; previous[v]=v; listdegree[d]=v;}
    else {next[v]=next[listdegree[d]]; next[listdegree[d]]=v;
      previous[v]=listdegree[d]; previous[next[v]]=v;} }

  vector<Node>removalorder(n);
  for (int i=0;i<n;i++) {// in each loop, we get one node of smallest degree
    Node v; int dv;
    for (dv=0;(dv<n)&&(listdegree[dv]==INVALID);dv++); // smallest degree is dv
    v = listdegree[dv];    removalorder[n-i-1] = v;    removed[v]=true;
    // remove the first node of degree dv
    if (listdegree[dv]==next[listdegree[dv]]){listdegree[dv]=INVALID;}//list with one node
    else {Node pre=previous[v], nex=next[v];
      next[pre] = nex;  previous[nex] = pre;
      if (listdegree[dv]==v) listdegree[dv]=nex;}
    // Now, decrease the degree of the adjacent nodes
    for (IncEdgeIt e(g,v); e!=INVALID; ++e) { Node u;
      if (g.u(e)==v) u = g.v(e); else u = g.u(e); // u is the adjacent node of v
      if (removed[u]) continue;
      int du=degree[u]; // remove u from the list listdegree[du]
      if (listdegree[du]==next[listdegree[du]]) {listdegree[du] = INVALID;}
      else {Node pre=previous[u], nex=next[u];
    next[pre] = nex;  previous[nex] = pre;
    if (listdegree[du]==u) listdegree[du]=nex;}
      degree[u]--;    du = degree[u];
      // insert u in the new degree position degree[u]
      if (listdegree[du]==INVALID){ next[u]=u; previous[u]=u; listdegree[du]=u;}
      else {next[u]=next[listdegree[du]]; next[listdegree[du]]=u;
    previous[u]=listdegree[du]; previous[next[u]]=u;}
    }} // At this point, we have the removal order in "removalorder"
  unsigned int maxremoveddegree=Max(g,degree); // maximum degree (considering degree when node was removed).
  int adjacentcolors[maxremoveddegree];
  for (NodeIt v(g);v!=INVALID;++v) nodecolor[v]=-1; // means no color
  for (int i=0;i<n;i++) {Node v =  removalorder[i]; // Let us paint node v
    int nadjcolors=0;
    for (IncEdgeIt e(g,v); e!=INVALID; ++e) { Node u;
      if (g.u(e)==v) u = g.v(e); else u = g.u(e); // u is the adjacent node of v
      if (nodecolor[u]!=-1) {adjacentcolors[nadjcolors]=nodecolor[u]; nadjcolors++;} }
    sort(adjacentcolors,adjacentcolors+nadjcolors);

    //cout<< "Cores adjacentes ao noh "<<vname[v]<<":  ";
    //for (int j=0;j<nadjcolors;j++) cout<<adjacentcolors[j]<<", ";   cout<<endl;

    { int j,color_v;
      if (nadjcolors==0) color_v=0;
      else if (adjacentcolors[0]!=0) color_v = 0;
      else if (nadjcolors==1)  color_v=1;
      else {// has at least two nodes and the first color is 0
    for (j=0;j<nadjcolors-1;j++) if (adjacentcolors[j+1]-adjacentcolors[j]>1) break;
    color_v=adjacentcolors[j]+1;}
      //cout << "Cor a ser atribuida ao vertice "<< vname[v]<<" "<<color_v<<endl<<endl;
      nodecolor[v] = color_v;
    }
  }
  ncolors = Max(g,nodecolor)+1;
  return(true);
}
