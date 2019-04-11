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
void readInstance(Digraph &dg, DNodePosMap &dposx, DNodePosMap &dposy, ArcValueMap &dweight, DNodeIntMap &disTerminal,
                  vector<DNode> &terminals, int &n, int &m, string filename);
#endif // MAIN_H
