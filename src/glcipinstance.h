#ifndef GLCIPINSTANCE_H
#define GLCIPINSTANCE_H

#include "mygraphlib.h"
#include <list>
#include <cassert>
#include "easyscip.h"
#include <utility>
#include <map>
#include <vector>

class GLCIPInstance
{
public:
    Digraph &g;
    DNodeStringMap &nodeName;
    DNodePosMap &posx;
    DNodePosMap &posy;
    ArcValueMap &influence;
    DNodeValueMap &threshold;
    DNodeValueVectorMap &incentives;
    double alpha;
    int n, m;

    GLCIPInstance(
        Digraph &_g, 
        DNodeStringMap &_nodeName, 
        DNodePosMap &_posx, 
        DNodePosMap &_posy, 
        ArcValueMap &_influence, 
        DNodeValueMap &_threshold,        
        DNodeValueVectorMap &_incentives, 
        double _alpha, 
        int _n, 
        int _m);
};

#endif // GLCIPINSTANCE_H
