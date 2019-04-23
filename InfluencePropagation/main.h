#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include "mygraphlib.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include "glcipinstance.h"
#include "mygraphlib.h"
#include "glcipsolution.h"
#include "arcmodel.h"
#include "graphviewer.h"

typedef struct structParams
{
    string alg;
    int timeLimit;
    string inputFile;
    bool graph;
    double alpha;
    //string    outputFile;
} Params;

void readCheckParams(Params &params, int argc, char *argv[]);
void readInstance(Digraph &g, DNodeStringMap &nodeName, DNodePosMap &posx, DNodePosMap &posy, ArcValueMap &influence, DNodeValueMap &threshold,
    DNodeValueVectorMap &incentives, int &n, int &m, string filename);
#endif // MAIN_H
