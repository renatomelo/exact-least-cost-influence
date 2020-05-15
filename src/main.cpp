#include "main.h"

int main(int argc, char *argv[])
{
    Digraph g;
    DNodeStringMap nodeName(g);
    DNodePosMap posx(g);
    DNodePosMap posy(g);
    ArcValueMap influence(g);
    DNodeValueMap threshold(g);
    DNodeValueVectorMap incentives(g);
    int n, m;

    // read parameters from command line
    Params params;
    readCheckParams(params, argc, argv);

    // read instance from file
    readInstance(g, nodeName, posx, posy, influence, threshold, incentives, n, m, "../in/" + params.inputFile);
    GLCIPInstance instance(g, nodeName, posx, posy, influence, threshold, incentives, params.alpha, n, m);

    // solve it
    GLCIPSolution solution(g);

    //auto started = chrono::high_resolution_clock::now();
    /*
    if(params.alg.compare("cov") == 0)
        GurobiBCP::run(instance, solution, measures, params.timeLimit);
    */
    if (params.alg.compare("arc") == 0)
    {
        ArcModel::run(instance, solution, params.timeLimit);
        //GraphViewer::ViewGLCIPSolution(instance, solution, "GLCIP Solution - Arc Model");
    }
    if (params.alg.compare("arcwb") == 0)
    {
        ArcModelWithBounds::run(instance, solution, params.timeLimit);
        //GraphViewer::ViewGLCIPSolution(instance, solution, "GLCIP Solution - Arc Model");
    }
    if (params.alg.compare("cov") == 0)
    {
        CovModelAllVariables::run(instance, solution, params.timeLimit);
        //GraphViewer::ViewGLCIPSolution(instance, solution, "GLCIP Solution - Cov Model with all variables");
    }
    if (params.alg.compare("covcg") == 0)
    {
        CovModel::run(instance, solution, params.timeLimit);
        //GraphViewer::ViewGLCIPSolution(instance, solution, "GLCIP Solution - Cov Model with column generation");
    }

    /* auto done = chrono::high_resolution_clock::now();
    int time = chrono::duration_cast<chrono::milliseconds>(done-started).count(); */

    return 0;
}

// read argv params
void readCheckParams(Params &params, int argc, char *argv[])
{
    params.alg = "";
    params.timeLimit = 1800;
    params.inputFile = "";
    params.graph = false;
    params.alpha = 0.5;

    // Read
    for (int i = 1; i < argc; i++)
    {
        const string arg(argv[i]);
        string next;

        if ((i + 1) < argc)
        {
            next = string(argv[i + 1]);
        }
        else
        {
            next = string("");
        }

        if (arg.find("-t") == 0 && next.size() > 0)
        {
            params.timeLimit = atoi(next.c_str());
            i++;
            continue;
        }

        if (arg.find("-i") == 0 && next.size() > 0)
        {
            params.inputFile = next;
            i++;
            continue;
        }

        if (arg.find("-alpha") == 0 && next.size() > 0)
        {
            params.alpha = atof(next.c_str());
            i++;
            continue;
        }

        if (arg.find("-a") == 0 && next.size() > 0)
        {
            params.alg = next;
            i++;
            continue;
        }

        if (arg.find("-g") == 0)
        {
            params.graph = true;
            continue;
        }

        cerr << "Invalid parameter." << endl;
        exit(1);
    }

    // Check
    if (params.inputFile == "")
    {
        cerr << "Input file should be specified" << endl;
        exit(1);
    }
}

// read file and create corresponding graph on the instance variable
void readInstance(Digraph &g,
                  DNodeStringMap &nodeName,
                  DNodePosMap &posx,
                  DNodePosMap &posy,
                  ArcValueMap &influence,
                  DNodeValueMap &threshold,
                  DNodeValueVectorMap &incentives,
                  int &n,
                  int &m,
                  string filename)
{

    // Read the graph (only nodes and edges)
    DigraphTable GT(filename, g);
    n = GT.getNDNodes();
    m = GT.getNArcs();

    // put zeros in the incentives
/*     for (DNodeIt v(g); v != INVALID; ++v)
    {
        incentives[v].push_back(0.0);
    } */

    // get nodes and arcs parameters
    GT.GetColumn("influence", influence);
    GT.GetColumn("threshold", threshold);
    //GT.GetColumn("incentives", incentives);
    GT.GetColumn("nodename", nodeName);
    GT.GetColumn("posx", posx);
    GT.GetColumn("posy", posy);

    //generate the incentives
    double maxThr = 0;
    for (DNodeIt v(g); v != INVALID; ++v)
    {
        //get max threshold
        maxThr = max(maxThr, threshold[v]);
    }

    for (DNodeIt v(g); v != INVALID; ++v)
    {
        incentives[v].push_back(0.0);
        incentives[v].push_back(0.25 * maxThr);
        incentives[v].push_back(0.5 * maxThr);
        incentives[v].push_back(0.75 * maxThr);
        incentives[v].push_back(maxThr);
    }
    
}
