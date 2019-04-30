// =============================================================
#ifndef MYGRAPHLIB_DEFINE
#define MYGRAPHLIB_DEFINE

#include<float.h>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<lemon/concepts/digraph.h>
#include<lemon/concepts/graph.h>
#include<lemon/concepts/maps.h>
#include<lemon/list_graph.h>
#include<lemon/gomory_hu.h>
#include<lemon/math.h>
#include<lemon/preflow.h>
#include<string>
#include "myutils.h"
#include "mycolor.h"
#include "geompack.h"
#include <limits.h>

using namespace std;
using namespace lemon;

class mygraphlibdefaultparameters {
    public:
        mygraphlibdefaultparameters();
        // Parameters used for graph visualization
        string graphvizdrawingprogram; // neato (default), twopi, circo, fdp, dot
        void setgraphviz(string prog);
        // Be careful to change this parameters, since the instances in the
        // folder ./instances/ are using these default terms for the header
        // Header table
        string inputgraphtable_nnodes; // "nnodes"
        string inputgraphtable_nedges; // "nedges", for the Graph type
        string inputgraphtable_narcs; // "narcs", for the Digraph type
        string inputgraphtable_type; // "type"
        // Node table
        string inputgraphtable_nodename; // "nodename"
        string inputgraphtable_posx; // "posx"
        string inputgraphtable_posy; // "posy"
        // Edge Table
        string inputgraphtable_endpoint1; // "endpoint1", for the Graph type
        string inputgraphtable_endpoint2; // "endpoint2", for the Graph type
        // Arc Table
        string inputgraphtable_tail; // "tail", for the Digraph type
        string inputgraphtable_head; // "head", for the Digraph type
} ;

inline mygraphlibdefaultparameters::mygraphlibdefaultparameters(){
    graphvizdrawingprogram = "neato";
    // Be careful to change this parameters, since the instances in the
    // folder ./instances/ are using these default terms to read graph files
    // Header table
    inputgraphtable_nnodes =  "nnodes";
    inputgraphtable_nedges = "nedges"; // for the Graph type
    inputgraphtable_narcs = "narcs";   // for the Digraph type
    inputgraphtable_type = "type";
    // Node table
    inputgraphtable_nodename = "nodename";
    inputgraphtable_posx = "posx";
    inputgraphtable_posy = "posy";
    // Edge Table
    inputgraphtable_endpoint1 = "endpoint1"; // for the Graph type
    inputgraphtable_endpoint2 = "endpoint2"; // for the Graph type
    // Arc Table
    inputgraphtable_tail = "tail"; // for the Digraph type
    inputgraphtable_head = "head"; // for the Digraph type
}

extern mygraphlibdefaultparameters mygraphlibparameters;

// Set Graphviz drawing program
inline void SetGraphviz(string prog){
    mygraphlibparameters.graphvizdrawingprogram = prog;
}

#define MAXLABELNAME 200
#define MAXLINE 1000
#define VIEW_DOT 0

typedef ListGraph Graph;
typedef Graph::Node Node;
typedef Graph::Edge Edge;

typedef Graph::NodeIt NodeIt;
typedef Graph::NodeMap<bool> CutMap;
typedef Graph::NodeMap<double> NodeValueMap;
typedef Graph::NodeMap<int> NodeColorMap;
typedef Graph::NodeMap<int> NodeIndexMap;
typedef Graph::NodeMap<int> NodeIntMap;
typedef Graph::NodeMap<bool> NodeBoolMap;
typedef Graph::NodeMap<string> NodeStringMap;
typedef Graph::NodeMap<double> NodePosMap;
typedef Graph::NodeMap<Node> NodeNodeMap;
typedef Graph::NodeMap<Edge> NodeEdgeMap;

typedef Graph::EdgeIt EdgeIt;
typedef Graph::IncEdgeIt IncEdgeIt;
typedef Graph::EdgeMap<double> EdgeValueMap;
typedef Graph::EdgeMap<int> EdgeColorMap;
typedef Graph::EdgeMap<int> EdgeIndexMap;
typedef Graph::EdgeMap<int> EdgeIntMap;
typedef Graph::EdgeMap<string> EdgeStringMap;
typedef Graph::EdgeMap<Edge> EdgeEdgeMap;
typedef Graph::EdgeMap<Node> EdgeNodeMap;
typedef Graph::EdgeMap<bool> EdgeBoolMap;
typedef Preflow<Graph, EdgeValueMap> PFType;

typedef ListDigraph Digraph;
typedef Digraph::Arc Arc;
typedef Digraph::ArcIt ArcIt;
typedef Digraph::InArcIt InArcIt;
typedef Digraph::OutArcIt OutArcIt;
typedef Digraph::ArcMap<double> ArcValueMap;
typedef Digraph::ArcMap<int> ArcIntMap;
typedef Digraph::ArcMap<string> ArcStringMap;
typedef Digraph::ArcMap<int> ArcColorMap;
typedef Digraph::ArcMap<bool> ArcBoolMap;

typedef Digraph::Node DNode;
typedef Digraph::NodeIt DNodeIt;
typedef Digraph::NodeMap<int> DNodeIntMap;
typedef Digraph::NodeMap<double> DNodeValueMap;
typedef Digraph::NodeMap<int> DNodeColorMap;
typedef Digraph::NodeMap<bool> DCutMap;
typedef Digraph::NodeMap<double> DNodePosMap;
typedef Digraph::NodeMap<vector<double>> DNodeValueVectorMap;
typedef Digraph::NodeMap<string> DNodeStringMap;
typedef Digraph::NodeMap<bool> DNodeBoolMap;
typedef Digraph::NodeMap<Arc> DNodeArcMap;
typedef Digraph::NodeMap<DNode> DNodeDNodeMap;

// used to storage a list of values for a node
typedef Digraph::NodeMap<vector<double>> DNodeListOfValuesMap;

// Used to obtain the node/dnode/edge/arc of a table
typedef vector<Node> LineToNodeMap;
typedef vector<DNode> LineToDNodeMap;
typedef vector<Edge> LineToEdgeMap;
typedef vector<Arc> LineToArcMap;

#if __cplusplus >= 201103L
#include <unordered_map>
typedef std::unordered_map<string,Node> StringToNodeMap;
typedef std::unordered_map<string,DNode> StringToDNodeMap;
#else
#include <tr1/unordered_map>
typedef std::tr1::unordered_map<string,Node> StringToNodeMap;
typedef std::tr1::unordered_map<string,DNode> StringToDNodeMap;
#endif


inline double Product(Graph &g,NodeValueMap &a,NodeValueMap &b)
{ double s=0.0; for (NodeIt v(g);v!=INVALID;++v) s += a[v]*b[v]; return(s);}
inline double Product(Graph &g,EdgeValueMap &a,EdgeValueMap &b)
{ double s=0.0; for (EdgeIt e(g);e!=INVALID;++e) s += a[e]*b[e]; return(s);}
inline double Product(Digraph &g,DNodeValueMap &a,DNodeValueMap &b)
{ double s=0.0; for (DNodeIt v(g);v!=INVALID;++v) s += a[v]*b[v]; return(s);}
inline double Product(Digraph &g,ArcValueMap &a,ArcValueMap &b)
{ double s=0.0; for (ArcIt e(g);e!=INVALID;++e) s += a[e]*b[e]; return(s);}


inline string GetGraphFileType(string filename){
  ifstream file;  file.open(filename.c_str());
  if (!file) {cout<<"Error: Could not open file "<<filename<<"."<<endl;exit(0);}
  StringTable Header(1,file);
  string type=Header.first("type");
  lowercase(type);
  file.close();
  return(type);}

bool ReadDigraph(string filename, Digraph &g, DNodeStringMap& vname,
         DNodePosMap& posx, DNodePosMap& posy, ArcValueMap& weight);
bool ReadGraph(string filename, Graph &g, NodeStringMap& vname,
           NodePosMap& posx, NodePosMap& posy, EdgeValueMap& weight);
bool ReadGraph(string filename, Graph &g, NodeStringMap& vname,
           NodePosMap& posx, NodePosMap& posy, NodeValueMap& weight);

bool InitGraphTable(StringTable &T,Graph &g,LineToNodeMap &nodemap);
bool InitGraphTable(StringTable &T,Graph &g,StringToNodeMap &smap,
            string colname_v1,string colname_v2,LineToEdgeMap &edgemap);
bool InitGraphTable(StringTable &T,Digraph &g,LineToDNodeMap &dnodemap);
bool InitGraphTable(StringTable &T,Digraph &g,StringToDNodeMap &smap,
            string colname_v1,string colname_v2,LineToArcMap &arcmap);

bool ReadGraphColumn(StringTable &T,LineToNodeMap &nodemap,string colname,NodeValueMap &read);
bool ReadGraphColumn(StringTable &T,LineToNodeMap &nodemap,string colname,NodeStringMap &read);
bool ReadGraphColumn(StringTable &T,LineToNodeMap &nodemap,string colname,NodeIntMap &read);
bool ReadGraphColumn(StringTable &T,LineToEdgeMap &edgemap,string colname,EdgeValueMap &read);
bool ReadGraphColumn(StringTable &T,LineToEdgeMap &edgemap,string colname,EdgeStringMap &read);
bool ReadGraphColumn(StringTable &T,LineToEdgeMap &edgemap,string colname,EdgeIntMap &read);
bool ReadGraphColumn(StringTable &T,LineToDNodeMap &nodemap,string colname,DNodeValueMap &read);
bool ReadGraphColumn(StringTable &T,LineToDNodeMap &nodemap,string colname,DNodeStringMap &read);
bool ReadGraphColumn(StringTable &T,LineToDNodeMap &nodemap,string colname,DNodeIntMap &read);
bool ReadGraphColumn(StringTable &T,LineToArcMap &arcmap,string colname,ArcValueMap &read);
bool ReadGraphColumn(StringTable &T,LineToArcMap &arcmap,string colname,ArcStringMap &read);
bool ReadGraphColumn(StringTable &T,LineToArcMap &arcmap,string colname,ArcIntMap &read);


// ==============================================================
bool GenerateVertexPositions(Graph   &g,NodePosMap  &posx, NodePosMap   &posy);
bool GenerateVertexPositions(Digraph &g,DNodePosMap &posx, DNodePosMap  &posy);
bool GenerateVertexPositions(Graph   &g,    EdgeValueMap &custo,
                 NodePosMap  &posx, NodePosMap   &posy);
bool GenerateVertexPositions(Digraph &g,    ArcValueMap  &custo,
                 DNodePosMap &posx, DNodePosMap  &posy);


void PrintGraph(Graph  &g, NodeStringMap &vname, EdgeValueMap &graphweight);



// Generate a triangulated Digraph, building the Delaunay
// triangulation of random points. Each edge of the Delaunay triangulation
// leads to two arcs (in both senses)
// Uses the geompack program, available in
// http://people.sc.fsu.edu/~jburkardt/cpp_src/geompack/geompack.html
bool GenerateTriangulatedDigraph(Digraph &g,
         DNodeStringMap &vname, // name of the nodes
         DNodePosMap    &px, // x-position of the nodes
         DNodePosMap    &py, // y-position of the nodes
         ArcValueMap     & weight, // weight of edges
         int n, // number of nodes
         double SizeX, // coordinate x is a random number in [0,SizeX)
                 double SizeY); // coordinate y is a random number in [0,SizeY)

// the same as above, but for non-oriented edges
bool GenerateTriangulatedGraph(Graph &g, // return with generated graph
              NodeStringMap &vname, // return with name of the nodes
              NodePosMap    & px, // return with x-position of the nodes
              NodePosMap    & py, // return with y-position of the nodes
              EdgeValueMap  & weight, // return with weight of edges
              int n, // number of nodes
              double SizeX, // coordinate x is a random number in [0,SizeX)
              double SizeY); // coordinate y is a random number in [0,SizeY)


class AdjacencyMatrix {
public:
  AdjacencyMatrix(Graph &graph,EdgeValueMap &graphweight,double NonEdgeValue);
  ~AdjacencyMatrix();
  double *AdjMatrix;
  Graph *g;
  EdgeValueMap *weight;
  int Nnodes,Nedges,Nmatrix;
  double NonEdgeValue;
  Node *Index2Node;
  Edge *Index2Edge;
  double Cost(Node,Node);
  double Cost(Edge);
  NodeIndexMap Node2Index;
  EdgeIndexMap Edge2Index;
};

//Generate a random complete euclidean Graph
bool GenerateRandomEuclideanGraph(Graph &g,
          NodeStringMap &vname, // node names
          NodePosMap    &px, // x-position of the nodes
          NodePosMap    &py, // y-position of the nodes
          EdgeValueMap  &weight, // weight of edges
          int n, // number of nodes
          double SizeX, // coordinate x is a random number in [0,SizeX)
              double SizeY); // coordinate y is a random number in [0,SizeY)

//Generate a random complete euclidean Digraph
bool GenerateRandomEuclideanDigraph(Digraph &g,
            DNodeStringMap &vname, // node names
            DNodePosMap    &px, // x-position of the nodes
            DNodePosMap    &py, // y-position of the nodes
            ArcValueMap  &weight, // weight of edges
            int n, // number of nodes
            double SizeX, // coordinate x is a random number in [0,SizeX)
            double SizeY); // coordinate y is a random number in [0,SizeY)


// Given a color code, return its name
//std::string ColorName(int cor);

double MinCut(Graph &g, EdgeValueMap &weight, Node &s,Node &t, CutMap &cut);


// Obtain a mininum cut for directed graphs from s to t.
// The returned cut is given by the vector of nodes 'cut' (boolean
// vector: nodes v in the same side of s have cut[v]=true, otherwise cut[v]=false.
double DiMinCut(Digraph &g,ArcValueMap &weight,DNode &s,DNode &t,DCutMap &cut);


Node GetNodeByName(Graph &g,
           NodeStringMap  &vname, // name of the nodes
           string vertexname);

// Return true if the edge vector is integer
inline bool IsInteger(Graph &g, EdgeValueMap &vx)
{for(EdgeIt e(g);e!=INVALID; ++e) if(IsFrac(vx[e]))return false; return(true);}


//--------------------------------------------------------------------------------
// The structure GraphTable is used to read a graph from a file.
// For Graph's, the file must have the following format,
// divided in three parts:
// #=========================================================
//    nnodes        nnedges         type          info1  info2 info3 ...
//  <num.nodes>  <num.edges>       graph           ... ... ...
// #---------------------------------------------------------
//  nodename nodeinfo1  nodeinfo2 ... ...
//  ...
//  ... (nnodes lines, one for each node)
//  ...
// #---------------------------------------------------------
//  endpoint1   endpoint2  edgeinfo1 edgeinfo2  edgeinfo3 ...
//  ...
//  ... (nnedges lines, one for each edge)
//  ...
// #=========================================================
//

// Required fields are:
//      nnodes, nnedges, type
//      nodename
//      endpoint1 endpoint2
// the other fields are optional. It depends on the informations you want to obtain.
// For some problems, there is only the edge weight.
// In other, there is node weight,...
class GraphTable {
    public:
        GraphTable(string filename,Graph &graph);
        ~GraphTable();
        int getNNodes();
        int getNEdges();
        bool GetColumn(string colname,NodeValueMap &col);
        bool GetColumn(string colname,NodeIntMap &col);
        bool GetColumn(string colname,NodeStringMap &col);
        bool GetColumn(string colname,EdgeValueMap &col);
        bool GetColumn(string colname,EdgeIntMap &col);
        bool GetColumn(string colname,EdgeStringMap &col);
        bool GetNodeCoordinates(string colx, NodePosMap &posx, string coly, NodePosMap &posy);
        void print();
        StringTable *Header,*NodeTable,*EdgeTable;
        Graph &g;
        LineToNodeMap line2node;
        LineToEdgeMap line2edge;
        StringToNodeMap name2node;
        int n;
        int m;
};

// Constructor to read graph given by list of nodes and
// list of edges
inline GraphTable::GraphTable(string filename,Graph &graph): g(graph){
    ifstream file;
    file.open(filename.c_str());

    NodeStringMap vname(this->g);

    if(!file){
      cout<<"Error: Could not open file "<<filename<<"."<<endl;
      exit(0);
    }

    this->Header = new StringTable(1,file);
    if (!this->Header){
      cout<<"Error: Memory allocation for Header problem.\n";
      exit(0);
    }

    string type = Header->first("type");
    lowercase(type);

    if(type == "graph"){
        int nnodes = StringToInt(Header->first("nnodes"));
        int nedges = StringToInt(Header->first("nedges"));
        this->n = nnodes;
        this->m = nedges;

        this->NodeTable = new StringTable(nnodes,file);
        if(!this->NodeTable){
            cout<<"Error: Memory allocation for NodeTable problem.\n";
            exit(0);
        }

        InitGraphTable(*(this->NodeTable), this->g, this->line2node);
        ReadGraphColumn(*(this->NodeTable),this->line2node,"nodename",vname);
        for (NodeIt v(g); v!=INVALID; ++v)
            name2node[vname[v]] = v;

        this->EdgeTable = new StringTable(nedges,file);
        if (!this->EdgeTable){
            cout<<"Error: Memory allocation for EdgeTable problem.\n";
            exit(0);
        }

        InitGraphTable(*(this->EdgeTable), this->g, this->name2node, "endpoint1", "endpoint2", this->line2edge);
    }

    file.close();
}

// Destructor
inline GraphTable::~GraphTable(){
    delete this->Header;
    delete this->NodeTable;
    delete this->EdgeTable;
}

inline bool GraphTable::GetColumn(string colname,NodeValueMap &column){
  int col = this->NodeTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->NodeTable->nrows;i++)
    this->NodeTable->entry(i,col,column[this->line2node[i]]);
  return(true);}
inline bool GraphTable::GetColumn(string colname,NodeIntMap &column){
  int col = this->NodeTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->NodeTable->nrows;i++)
    this->NodeTable->entry(i,col,column[this->line2node[i]]);
  return(true);}
inline bool GraphTable::GetColumn(string colname,NodeStringMap &column){
  int col = this->NodeTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->NodeTable->nrows;i++)
    this->NodeTable->entry(i,col,column[this->line2node[i]]);
  return(true);}
inline bool GraphTable::GetColumn(string colname,EdgeValueMap &column){
  int col = this->EdgeTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->EdgeTable->nrows;i++)
    this->EdgeTable->entry(i,col,column[this->line2edge[i]]);
  return(true);}
inline bool GraphTable::GetColumn(string colname,EdgeIntMap &column){
  int col = this->EdgeTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->EdgeTable->nrows;i++)
    this->EdgeTable->entry(i,col,column[this->line2edge[i]]);
  return(true);}
inline bool GraphTable::GetColumn(string colname,EdgeStringMap &column){
  int col = this->EdgeTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->EdgeTable->nrows;i++)
    this->EdgeTable->entry(i,col,column[this->line2edge[i]]);
  return(true);}

inline void GraphTable::print(){
  this->Header->print();
  this->NodeTable->print();
  this->EdgeTable->print();
}

inline int GraphTable::getNNodes(){
    return this->n;
}

inline int GraphTable::getNEdges(){
    return this->m;
}

//--------------------------------------------------------------------------------
// This structure is used to read a digraph from a file.
// For Digraph's, the file must have the following format,
// divided in three parts:
// #=========================================================
//    nnodes        nnarcs         type          info1  info2 info3 ...
//  <num.nodes>  <num.edges>      digraph  ... ... ...
// #---------------------------------------------------------
//  nodename nodeinfo1  nodeinfo2 ... ...
//  ...
//  ... (nnodes lines, one for each node)
//  ...
// #---------------------------------------------------------
//  tail   head  arcinfo1 arcinfo2  arcinfo3 ...
//  ...
//  ... (nnarcs lines, one for each arc)
//  ...
// #=========================================================
//

// Required fields are:
//      nnodes, nnedges, type
//      nodename
//      endpoint1 endpoint2
// the other fields are optional. It depends the informations you want to obtain.
// For some problems, there is only the edge weight.
// In other, there is node weight,...
class DigraphTable {
public:
  DigraphTable(string filename,Digraph &graph);
  ~DigraphTable();
  int getNDNodes();
  int getNArcs();
  bool GetColumn(string colname,DNodeValueMap &col);
  bool GetColumn(string colname,DNodeIntMap &col);
  bool GetColumn(string colname,DNodeStringMap &col);
  bool GetColumn(string colname,ArcValueMap &col);
  bool GetColumn(string colname,ArcIntMap &col);
  bool GetColumn(string colname,ArcStringMap &col);
  bool GetColumn(string colname, DNodeValueVectorMap &column);
  void print();
  StringTable *Header,*NodeTable,*ArcTable;
  Digraph &g;
  LineToDNodeMap line2node;
  LineToArcMap line2arc;
  StringToDNodeMap name2node;
  int n;
  int m;
};
// Constructor to read graph given by list of nodes and
// list of arcs
inline DigraphTable::DigraphTable(string filename,Digraph &graph):
  g(graph)
{
  ifstream file;  file.open(filename.c_str());
  DNodeStringMap vname(this->g);
  if (!file) {cout<<"Error: Could not open file "<<filename<<"."<<endl;exit(0);}

  this->Header = new StringTable(1,file);
  if (!this->Header) {cout<<"Error: Memory allocation for Header problem.\n";exit(0);}
  string type=Header->first("type");
  lowercase(type);
  if(type=="digraph"){
    int nnodes=StringToInt(Header->first("nnodes"));
    int narcs=StringToInt(Header->first("narcs"));
    this->n = nnodes;
    this->m = narcs;

    this->NodeTable = new StringTable(nnodes,file);
    if (!this->NodeTable) {cout<<"Error: Memory allocation for NodeTable problem.\n";exit(0);}
    InitGraphTable(*(this->NodeTable),this->g,this->line2node);
    ReadGraphColumn(*(this->NodeTable),this->line2node,"nodename",vname);
    for (DNodeIt v(g); v!=INVALID; ++v) name2node[vname[v]] = v;

    this->ArcTable = new StringTable(narcs,file);
    if (!this->ArcTable) {cout<<"Error: Memory allocation for ArcTable problem.\n";exit(0);}
    InitGraphTable(*(this->ArcTable),this->g,this->name2node,"tail","head",this->line2arc);
  }
  file.close();
}

// Destructor
inline DigraphTable::~DigraphTable()
{ delete this->Header;  delete this->NodeTable; delete this->ArcTable;}


inline bool DigraphTable::GetColumn(string colname,DNodeValueMap &column){
  int col = this->NodeTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->NodeTable->nrows;i++)
    this->NodeTable->entry(i,col,column[this->line2node[i]]);
  return(true);}
inline bool DigraphTable::GetColumn(string colname,DNodeIntMap &column){
  int col = this->NodeTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->NodeTable->nrows;i++)
    this->NodeTable->entry(i,col,column[this->line2node[i]]);
  return(true);}
inline bool DigraphTable::GetColumn(string colname,DNodeStringMap &column){
  int col = this->NodeTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->NodeTable->nrows;i++)
    this->NodeTable->entry(i,col,column[this->line2node[i]]);
  return(true);}
inline bool DigraphTable::GetColumn(string colname,ArcValueMap &column){
  int col = this->ArcTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->ArcTable->nrows;i++)
    this->ArcTable->entry(i,col,column[this->line2arc[i]]);
  return(true);}
inline bool DigraphTable::GetColumn(string colname,ArcIntMap &column){
  int col = this->ArcTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->ArcTable->nrows;i++)
    this->ArcTable->entry(i,col,column[this->line2arc[i]]);
  return(true);}
inline bool DigraphTable::GetColumn(string colname,ArcStringMap &column){
  int col = this->ArcTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->ArcTable->nrows;i++)
    this->ArcTable->entry(i,col,column[this->line2arc[i]]);
  return(true);}
inline bool DigraphTable::GetColumn(string colname, DNodeValueVectorMap &column){
  int col = this->NodeTable->columnindex(colname);
  if (col==-1) {cout << "Error: Column name \""+colname+"\" missing.\n";exit(0);}
  for (int i=0;i<this->NodeTable->nrows;i++){
    string entry;
    this->NodeTable->entry(i,col,entry);

    std::string delim = ",";
    auto start = 0;
    auto end = entry.find(delim);
    while (end != std::string::npos)
    {
      column[this->line2node[i]].push_back(stod(entry.substr(start, end - start)));
      start = end + delim.length();
      end = entry.find(delim, start);
    }

    column[this->line2node[i]].push_back(stod(entry.substr(start, end)));
  }
  return(true);}

inline void DigraphTable::print(){
  this->Header->print();
  this->NodeTable->print();
  this->ArcTable->print();
}

inline int DigraphTable::getNDNodes(){
    return this->n;
}

inline int DigraphTable::getNArcs(){
    return this->m;
}

bool GetNodeCoordinates(DigraphTable &GT, string colx, DNodePosMap &posx, string coly, DNodePosMap &posy);

//--------------------------------------------------------------------------------
class GraphAttributes {
public:
  GraphAttributes(Graph &graph, NodeStringMap &vname,
          NodePosMap &posx, NodePosMap &posy);
  GraphAttributes(GraphTable &GT);

  Graph &g;

  // Graph attributes
  string    graph_attrib; // Graph attributes
  string    default_node_attrib; //  Default node attributes
  string    default_edge_attrib; //  Default edge attributes

  NodeStringMap  node_name;
  NodePosMap     px;
  NodePosMap     py;
  NodeStringMap  node_attrib;
  EdgeStringMap  edge_attrib;
  double factor;
  void SetImageFactor(double f); // 1.0 does not scale the objects
  void SetGraphAttrib(string s);
  void SetDefaultNodeAttrib(string s);
  void SetDefaultEdgeAttrib(string s);
  void SetAttrib(Node v,string s);
  void SetAttrib(Edge e,string s);
  void SetShape(Node v,string s);
  void SetShape(Edge e,string s);
  void SetStyle(Node v,string s);
  void SetStyle(Edge e,string s);
  void SetColor(Node v,string s);
  void SetColor(Edge e,string s);
  void SetLabel(Node v,string s);
  void SetLabel(Edge e,string s);
  void SetLabel(Node v,double d);
  void SetLabel(Edge e,double d);
  void SetLabel(string s);
  void SetLabel(EdgeValueMap &x);
  void SetLabelNonZero(EdgeValueMap &x);
  void SetLabel(NodeValueMap &x);
  void SetLabelNonZero(NodeValueMap &x);
  void SetLabel(EdgeStringMap &x);
  void SetLabel(NodeStringMap &x);

  void SetDefaultNodeFontName(string s);
  void SetDefaultEdgeFontName(string s);
  void SetFontName(Node v, string s);
  void SetFontName(Edge e, string s);
  void SetDefaultNodeFontSize(double s);
  void SetDefaultEdgeFontSize(double s);
  void SetFontSize(Node v, double s);
  void SetFontSize(Edge e, double s);
  void SetFontSize(double s);

  void SetColorByValue(int value, NodeIntMap &nodeval, string cor);
  void SetColorByValue(double value, NodeValueMap &nodeval, string cor);
  void SetColorByInterval(double lb,double ub, NodeValueMap &nodeval, string cor);
  void SetColorByValue(double value, EdgeValueMap &edgeval, string cor);
  void SetColorByInterval(double lb,double ub, EdgeValueMap &edgeval, string cor);

  void SetAttribByValue(int value, NodeIntMap &nodeval, string at);
  void SetAttribByValue(double value, NodeValueMap &nodeval, string at);
  void SetAttribByInterval(double lb,double ub, NodeValueMap &nodeval, string at);
  void SetAttribByValue(int value, EdgeIntMap &edgeval, string at);
  void SetAttribByValue(double value, EdgeValueMap &edgeval, string at);
  void SetAttribByInterval(double lb,double ub, EdgeValueMap &edgeval, string at);

  bool Save(string name);
  bool View();
};


inline void GraphAttributes::SetImageFactor(double f){ this->factor = f;}

inline void GraphAttributes::SetGraphAttrib(string s){ this->graph_attrib += " "+s;}

inline void GraphAttributes::SetDefaultNodeAttrib(string s){ this->default_node_attrib += " "+s; }
inline void GraphAttributes::SetDefaultEdgeAttrib(string s){ this->default_edge_attrib += " "+s;}
inline void GraphAttributes::SetAttrib(Node v,string s){ this->node_attrib[v] += " "+s;}
inline void GraphAttributes::SetAttrib(Edge e,string s){ this->edge_attrib[e] += " "+s;}

inline void GraphAttributes::SetShape(Node v,string s){this->SetAttrib(v,"shape=\""+s+"\"");}
inline void GraphAttributes::SetShape(Edge e,string s){this->SetAttrib(e,"shape=\""+s+"\"");}

inline void GraphAttributes::SetStyle(Node v,string s){this->SetAttrib(v,"style=\""+s+"\"");}
inline void GraphAttributes::SetStyle(Edge e,string s){this->SetAttrib(e,"style=\""+s+"\"");}

inline void GraphAttributes::SetColor(Node v,string s){
  this->SetAttrib(v,"color=\""+s+"\"");}

inline void GraphAttributes::SetColor(Edge e,string s){
  this->SetAttrib(e,"color=\""+s+"\"");}


inline void GraphAttributes::SetLabel(Node v,string s){this->SetAttrib(v,"label=\""+s+"\"");}
inline void GraphAttributes::SetLabel(Edge e,string s){this->SetAttrib(e,"label=\""+s+"\"");}
inline void GraphAttributes::SetLabel(Node v,double d){this->SetLabel(v,DoubleToString(d));}
inline void GraphAttributes::SetLabel(Edge e,double d){this->SetLabel(e,DoubleToString(d));}
inline void GraphAttributes::SetLabel(string s){this->SetGraphAttrib("label=\""+s+"\"");}

inline void GraphAttributes::SetLabel(EdgeValueMap &x)
{ for (EdgeIt e(this->g);e!=INVALID;++e) this->SetLabel(e,x[e]);}
inline void GraphAttributes::SetLabelNonZero(EdgeValueMap &x)
{ for (EdgeIt e(this->g);e!=INVALID;++e) if(!IsEqual(0,x[e])) this->SetLabel(e,x[e]);}
inline void GraphAttributes::SetLabel(NodeValueMap &x)
{ for (NodeIt v(this->g);v!=INVALID;++v) this->SetLabel(v,x[v]);}
inline void GraphAttributes::SetLabelNonZero(NodeValueMap &x)
{ for (NodeIt v(this->g);v!=INVALID;++v) if(!IsEqual(0,x[v])) this->SetLabel(v,x[v]);}
inline void GraphAttributes::SetLabel(EdgeStringMap &x)
{ for (EdgeIt a(this->g);a!=INVALID;++a) this->SetLabel(a,x[a]);}
inline void GraphAttributes::SetLabel(NodeStringMap &x)
{ for (NodeIt v(this->g);v!=INVALID;++v) this->SetLabel(v,x[v]);}


inline void GraphAttributes::SetDefaultNodeFontName(string s){this->SetDefaultNodeAttrib("fontname=\""+s+"\"");}
inline void GraphAttributes::SetDefaultEdgeFontName(string s){this->SetDefaultEdgeAttrib("fontname=\""+s+"\"");}
inline void GraphAttributes::SetFontName(Node v, string s){this->SetAttrib(v,"fontname=\""+s+"\"");}
inline void GraphAttributes::SetFontName(Edge e, string s){this->SetAttrib(e,"fontname=\""+s+"\"");}

inline void GraphAttributes::SetDefaultNodeFontSize(double s){this->SetDefaultNodeAttrib("fontsize=\""+DoubleToString(s)+"\"");}
inline void GraphAttributes::SetDefaultEdgeFontSize(double s){this->SetDefaultEdgeAttrib("fontsize=\""+DoubleToString(s)+"\"");}
inline void GraphAttributes::SetFontSize(Node v, double s){this->SetAttrib(v,"fontsize=\""+DoubleToString(s)+"\"");}
inline void GraphAttributes::SetFontSize(Edge e, double s){this->SetAttrib(e,"fontsize=\""+DoubleToString(s)+"\"");}
inline void GraphAttributes::SetFontSize(double s){this->SetGraphAttrib("fontsize=\""+DoubleToString(s)+"\"");}


inline void GraphAttributes::SetColorByValue(int value,
                      NodeIntMap &nodeval, string cor){

  for (NodeIt v(this->g);v!=INVALID;++v) if(value==nodeval[v]) this->SetColor(v,cor);}

inline void GraphAttributes::SetColorByValue(double value,
                      NodeValueMap &nodeval, string cor){

  for (NodeIt v(this->g);v!=INVALID;++v) if(IsEqual(value,nodeval[v])) this->SetColor(v,cor);}

inline void GraphAttributes::SetColorByInterval(double lb,double ub,
                      NodeValueMap &nodeval, string cor){
  for (NodeIt v(this->g);v!=INVALID;++v)
    if((nodeval[v]>lb-MY_EPS) && (nodeval[v]<ub+MY_EPS)) this->SetColor(v,cor);}

inline void GraphAttributes::SetColorByValue(double value,
                      EdgeValueMap &edgeval, string cor){
  for (EdgeIt e(this->g);e!=INVALID;++e) if(IsEqual(value,edgeval[e])) this->SetColor(e,cor);}

inline void GraphAttributes::SetColorByInterval(double lb,double ub,
                      EdgeValueMap &edgeval, string cor){
  for (EdgeIt e(this->g);e!=INVALID;++e)
    if((edgeval[e]>lb-MY_EPS) && (edgeval[e]<ub+MY_EPS)) this->SetColor(e,cor);}


inline void GraphAttributes::SetAttribByValue(int value,
                      NodeIntMap &nodeval, string at){

  for (NodeIt v(this->g);v!=INVALID;++v) if(value==nodeval[v]) this->SetAttrib(v,at);}

inline void GraphAttributes::SetAttribByValue(double value,
                      NodeValueMap &nodeval, string at){

  for (NodeIt v(this->g);v!=INVALID;++v) if(IsEqual(value,nodeval[v])) this->SetAttrib(v,at);}

inline void GraphAttributes::SetAttribByInterval(double lb,double ub,
                      NodeValueMap &nodeval, string at){
  for (NodeIt v(this->g);v!=INVALID;++v)
    if((nodeval[v]>lb-MY_EPS) && (nodeval[v]<ub+MY_EPS)) this->SetAttrib(v,at);}

inline void GraphAttributes::SetAttribByValue(int value,
                      EdgeIntMap &edgeval, string at){
  for (EdgeIt e(this->g);e!=INVALID;++e) if(value==edgeval[e]) this->SetAttrib(e,at);}

inline void GraphAttributes::SetAttribByValue(double value,
                      EdgeValueMap &edgeval, string at){
  for (EdgeIt e(this->g);e!=INVALID;++e) if(IsEqual(value,edgeval[e])) this->SetAttrib(e,at);}

inline void GraphAttributes::SetAttribByInterval(double lb,double ub,
                      EdgeValueMap &edgeval, string at){
  for (EdgeIt e(this->g);e!=INVALID;++e)
    if((edgeval[e]>lb-MY_EPS) && (edgeval[e]<ub+MY_EPS)) this->SetAttrib(e,at);}




//--------------------------------------------------------------------------------
class DigraphAttributes {
public:
  DigraphAttributes(Digraph &digraph, DNodeStringMap &vname,
          DNodePosMap &posx, DNodePosMap &posy);
  Digraph &g;

  // Graph attributes
  string    digraph_attrib; // Graph attributes
  string    default_node_attrib; //  Default node attributes
  string    default_arc_attrib; //  Default arc attributes

  DNodeStringMap  &node_name;
  DNodePosMap     &px;
  DNodePosMap     &py;
  DNodeStringMap  node_attrib;
  ArcStringMap  arc_attrib;
  double factor;
  void SetImageFactor(double f); // 1.0 does not scale the objects
  void SetDigraphAttrib(string s);
  void SetDefaultDNodeAttrib(string s);
  void SetDefaultArcAttrib(string s);
  void SetAttrib(DNode v,string s);
  void SetAttrib(Arc e,string s);

  void SetShape(DNode v,string s);
  void SetShape(Arc e,string s);
  void SetStyle(DNode v,string s);
  void SetStyle(Arc e,string s);
  void SetColor(DNode v,string s);
  void SetColor(Arc e,string s);
  void SetLabel(DNode v,string s);
  void SetLabel(Arc e,string s);
  void SetLabel(DNode v,double d);
  void SetLabel(Arc e,double d);
  void SetLabel(string s);
  void SetLabel(ArcValueMap &x);
  void SetLabel(DNodeValueMap &x);
  void SetLabel(ArcStringMap &x);
  void SetLabel(DNodeStringMap &x);

  void SetDefaultDNodeFontName(string s);
  void SetDefaultArcFontName(string s);
  void SetFontName(DNode v, string s);
  void SetFontName(Arc e, string s);

  void SetDefaultDNodeFontSize(double s);
  void SetDefaultArcFontSize(double s);
  void SetFontSize(DNode v, double s);
  void SetFontSize(Arc e, double s);
  void SetFontSize(double s);

  void SetColorByValue(double value, DNodeValueMap &nodeval, string cor);
  void SetColorByValue(double value, ArcValueMap &arcval, string cor);
  void SetColorByInterval(double lb,double ub,DNodeValueMap &nodeval, string cor);
  void SetColorByInterval(double lb,double ub,ArcValueMap &arcval, string cor);

  void SetAttribByValue(int value, DNodeIntMap &nodeval, string at);
  void SetAttribByValue(double value, DNodeValueMap &nodeval, string at);
  void SetAttribByValue(int value, ArcIntMap &arcval, string at);
  void SetAttribByValue(double value, ArcValueMap &arcval, string at);
  void SetAttribByInterval(double lb,double ub,DNodeValueMap &nodeval, string at);
  void SetAttribByInterval(double lb,double ub,ArcValueMap &arcval, string at);


  bool View(); // text displayed below the figure
};
inline void DigraphAttributes::SetImageFactor(double f){ this->factor = f;}

inline void DigraphAttributes::SetDigraphAttrib(string s){ this->digraph_attrib += " "+s;}

inline void DigraphAttributes::SetDefaultDNodeAttrib(string s){ this->default_node_attrib += " "+s; }
inline void DigraphAttributes::SetDefaultArcAttrib(string s){ this->default_arc_attrib += " "+s;}
inline void DigraphAttributes::SetAttrib(DNode v,string s){ this->node_attrib[v] += " "+s;}
inline void DigraphAttributes::SetAttrib(Arc e,string s){ this->arc_attrib[e] += " "+s;}

inline void DigraphAttributes::SetShape(DNode v,string s){this->SetAttrib(v,"shape=\""+s+"\"");}
inline void DigraphAttributes::SetShape(Arc e,string s){this->SetAttrib(e,"shape=\""+s+"\"");}

inline void DigraphAttributes::SetStyle(DNode v,string s){this->SetAttrib(v,"style=\""+s+"\"");}
inline void DigraphAttributes::SetStyle(Arc e,string s){this->SetAttrib(e,"style=\""+s+"\"");}

inline void DigraphAttributes::SetColor(DNode v,string s){
  this->SetAttrib(v,"color=\""+s+"\"");}

inline void DigraphAttributes::SetColor(Arc e,string s){
  this->SetAttrib(e,"color=\""+s+"\"");}

inline void DigraphAttributes::SetLabel(DNode v,string s){this->SetAttrib(v,"label=\""+s+"\"");}
inline void DigraphAttributes::SetLabel(Arc e,string s){this->SetAttrib(e,"label=\""+s+"\"");}
inline void DigraphAttributes::SetLabel(DNode v,double d){this->SetLabel(v,DoubleToString(d));}
inline void DigraphAttributes::SetLabel(Arc e,double d){this->SetLabel(e,DoubleToString(d));}
inline void DigraphAttributes::SetLabel(string s){this->SetDigraphAttrib("label=\""+s+"\"");}
inline void DigraphAttributes::SetLabel(ArcValueMap &x)
{ for (ArcIt a(this->g);a!=INVALID;++a) this->SetLabel(a,x[a]);}
inline void DigraphAttributes::SetLabel(DNodeValueMap &x)
{ for (DNodeIt v(this->g);v!=INVALID;++v) this->SetLabel(v,x[v]);}
inline void DigraphAttributes::SetLabel(ArcStringMap &x)
{ for (ArcIt a(this->g);a!=INVALID;++a) this->SetLabel(a,x[a]);}
inline void DigraphAttributes::SetLabel(DNodeStringMap &x)
{ for (DNodeIt v(this->g);v!=INVALID;++v) this->SetLabel(v,x[v]);}



inline void DigraphAttributes::SetDefaultDNodeFontName(string s){this->SetDefaultDNodeAttrib("fontname=\""+s+"\"");}
inline void DigraphAttributes::SetDefaultArcFontName(string s){this->SetDefaultArcAttrib("fontname=\""+s+"\"");}
inline void DigraphAttributes::SetFontName(DNode v, string s){this->SetAttrib(v,"fontname=\""+s+"\"");}
inline void DigraphAttributes::SetFontName(Arc e, string s){this->SetAttrib(e,"fontname=\""+s+"\"");}

inline void DigraphAttributes::SetDefaultDNodeFontSize(double s){this->SetDefaultDNodeAttrib("fontsize=\""+DoubleToString(s)+"\"");}
inline void DigraphAttributes::SetDefaultArcFontSize(double s ){this->SetDefaultArcAttrib("fontsize=\""+DoubleToString(s)+"\"");}
inline void DigraphAttributes::SetFontSize(DNode v, double s){this->SetAttrib(v,"fontsize=\""+DoubleToString(s)+"\"");}
inline void DigraphAttributes::SetFontSize(Arc e, double s){this->SetAttrib(e,"fontsize=\""+DoubleToString(s)+"\"");}
inline void DigraphAttributes::SetFontSize(double s){this->SetDigraphAttrib("fontsize=\""+DoubleToString(s)+"\"");}

// ----
inline void DigraphAttributes::SetColorByValue(double value,DNodeValueMap &nodeval, string cor){
  for (DNodeIt v(this->g);v!=INVALID;++v) if(IsEqual(value,nodeval[v])) this->SetColor(v,cor);}

inline void DigraphAttributes::SetColorByInterval(double lb,double ub,DNodeValueMap &nodeval, string cor){
  for (DNodeIt v(this->g);v!=INVALID;++v)
    if((nodeval[v]>lb-MY_EPS) && (nodeval[v]<ub+MY_EPS)) this->SetColor(v,cor);}

inline void DigraphAttributes::SetColorByValue(double value, ArcValueMap &arcval, string cor){
  for (ArcIt e(this->g); e!=INVALID; ++e) if(IsEqual(value,arcval[e])) this->SetColor(e,cor);}

inline void DigraphAttributes::SetColorByInterval(double lb,double ub,ArcValueMap &arcval, string cor){
  for (ArcIt e(this->g); e!=INVALID; ++e) if((arcval[e]>lb-MY_EPS) && (arcval[e]<ub+MY_EPS)) this->SetColor(e,cor);}
// ----
inline void DigraphAttributes::SetAttribByValue(double value,DNodeValueMap &nodeval, string at){
  for (DNodeIt v(this->g);v!=INVALID;++v) if(IsEqual(value,nodeval[v])) this->SetAttrib(v,at);}

inline void DigraphAttributes::SetAttribByValue(int value,DNodeIntMap &nodeval, string at){
  for (DNodeIt v(this->g);v!=INVALID;++v) if(value==nodeval[v]) this->SetAttrib(v,at);}

inline void DigraphAttributes::SetAttribByInterval(double lb,double ub,DNodeValueMap &nodeval, string at){
  for (DNodeIt v(this->g);v!=INVALID;++v)
    if((nodeval[v]>lb-MY_EPS) && (nodeval[v]<ub+MY_EPS)) this->SetAttrib(v,at);}

inline void DigraphAttributes::SetAttribByValue(double value, ArcValueMap &arcval, string at){
  for (ArcIt e(this->g); e!=INVALID; ++e) if(IsEqual(value,arcval[e])) this->SetAttrib(e,at);}

inline void DigraphAttributes::SetAttribByValue(int value, ArcIntMap &arcval, string at){
  for (ArcIt e(this->g); e!=INVALID; ++e) if(value==arcval[e]) this->SetAttrib(e,at);}

inline void DigraphAttributes::SetAttribByInterval(double lb,double ub,ArcValueMap &arcval, string at){
  for (ArcIt e(this->g); e!=INVALID; ++e) if((arcval[e]>lb-MY_EPS) && (arcval[e]<ub+MY_EPS)) this->SetAttrib(e,at);}


// ---

typedef Graph::NodeMap<DNode> NodeDNodeMap;
typedef Graph::EdgeMap<Arc> EdgeArcMap;
bool Graph2Digraph(Graph &g, Digraph &dg, NodeDNodeMap &nodemap, EdgeArcMap &edgemap1, EdgeArcMap &edgemap2);

//bool greedycoloring(Graph &g,int &ncolors,NodeIntMap &nodecolor);
bool GreedyColoring(Graph &g,int &ncolors,NodeIntMap &nodecolor);


//Maximum and Minimum values of a vector
inline double Max(Graph &g,NodeValueMap &val){
  double M=-DBL_MAX; for (NodeIt v(g);v!=INVALID;++v) if (val[v]>M) M=val[v]; return(M);}
inline int Max(Graph &g,NodeIntMap &val){
  int M=INT_MIN; for (NodeIt v(g);v!=INVALID;++v) if (val[v]>M) M=val[v]; return(M);}
inline double Min(Graph &g,NodeValueMap &val){
  double M=DBL_MAX; for (NodeIt v(g);v!=INVALID;++v) if (val[v]<M) M=val[v]; return(M);}
inline int Min(Graph &g,NodeIntMap &val){
  int M=INT_MAX; for (NodeIt v(g);v!=INVALID;++v) if (val[v]<M) M=val[v]; return(M);}

#endif
