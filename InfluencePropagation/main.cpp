#include "main.h"

int main(int argc, char *argv[]){
    Digraph g;
    DNodePosMap posx(g);
    DNodePosMap posy(g);
    ArcValueMap influence;
    DNodeValueMap threshold;
    DNodeValueMap incentives;
    DNodeValueMap costs;
    int n, m;

    // read parameters from command line
    Params params;
    readCheckParams(params, argc, argv);

    // read instance from file
    readInstance(g, posx, posy, influence, threshold, incentives, costs, n, m, "in/" + params.inputFile);
    GLCIPInstance instance(g, posx, posy, influence, threshold, incentives, costs, params.alpha, n, m);

    // solve it
    GLCIPSolution solution(g);

    auto started = chrono::high_resolution_clock::now();
    Measures measures;
    if(params.alg.compare("cov") == 0)
        GurobiBCP::run(instance, solution, measures, params.timeLimit);
    if(params.alg.compare("arc") == 0)
        SCIPBCP::run(instance, solution, measures, params.timeLimit);

    auto done = chrono::high_resolution_clock::now();
    int time = chrono::duration_cast<chrono::milliseconds>(done-started).count();

    return 0;
}

// read argv params
void readCheckParams(Params &params,int argc, char *argv[])
{
    params.alg        = "";
    params.timeLimit  = 3600;
    params.inputFile  = "";
    params.graph = false;
    params.initialSolution = false;
    params.shouldPrice = false;
    params.alpha = 0.5;

    // Read
    for(int i = 1; i < argc; i++){
        const string arg(argv[i]);
        string next;

        if((i+1) < argc){
            next = string(argv[i+1]);
        }
        else{
            next = string("");
        }

        if(arg.find("-t") == 0 && next.size() > 0){
            params.timeLimit = atoi(next.c_str());
            i++;
            continue;
        }

        if(arg.find("-i") == 0 && next.size() > 0){
            params.inputFile = next;
            i++;
            continue;
        }

        if(arg.find("-a") == 0 && next.size() > 0){
            params.alg = next;
            i++;
            continue;
        }

        if(arg.find("-g") == 0){
            params.graph = true;
            continue;
        }

        if(arg.find("-alpha") == 0){
            params.alpha = atof(next.c_str());
            continue;
        }

        cerr << "Invalid parameter." << endl;
        exit(1);
    }

    // Check
    if(params.inputFile == ""){
        cerr << "Input file should be specified" << endl;
        exit(1);
    }
}

// read file and create corresponding graph on the instance variable
void readInstance(Digraph &g, DNodePosMap &posx, DNodePosMap &posy, ArcValueMap &influence, DNodeValueMap &threshold,
    DNodeValueMap &incentives, DNodeValueMap &costs, int &n, int &m, string filename){

    // Read the graph (only nodes and edges)
    GraphTable GT(filename, g);
    n = GT.getNNodes();
    m = GT.getNEdges();

    // get nodes and arcs parameters
    GT.GetColumn("influece", influence);
    GT.GetColumn("threshold", threshold);
    GT.GetColumn("incentives", incentives);
    GT.GetColumn("costs", costs);
    GT.GetNodeCoordinates("posx", posx, "posy", posy);
}
