#ifndef GLCIPINSTANCE_H
#define GLCIPINSTANCE_H

#include "mygraphlib.h"
#include <list>
#include <cassert>
#include "easyscip.h"
#include <utility>
#include <map>
#include <vector>
#include "star.h"

class GLCIPInstance
{
public:
    Digraph &g;
    DNodePosMap &posx;
    DNodePosMap &posy;
    ArcValueMap &influence;
    DNodeValueMap &threshold;
    DNodeValueVectorMap &incentives;
    double alpha;
    int n, m;

    GLCIPInstance(Digraph &g, DNodePosMap &posx, DNodePosMap &posy, ArcValueMap &influence, DNodeValueMap &threshold,
        DNodeValueVectorMap &incentives, double alpha, int n, int m);
};

#endif // GLCIPINSTANCE_H
