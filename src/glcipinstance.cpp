#include "glcipinstance.h"

GLCIPInstance::GLCIPInstance(
    Digraph &_g, 
    DNodeStringMap &_nodeName, 
    DNodePosMap &_posx, 
    DNodePosMap &_posy, 
    ArcValueMap &_influence,
    DNodeValueMap &_threshold, 
    DNodeValueVectorMap &_incentives, 
    double _alpha, 
    int _n, 
    int _m) : g(_g), 
    nodeName(_nodeName), 
    posx(_posx), 
    posy(_posy),
    influence(_influence),
    threshold(_threshold), 
    incentives(_incentives)
{
    this->alpha = _alpha;
    this->n = _n;
    this->m = _m;
}
