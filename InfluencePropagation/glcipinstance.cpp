#include "glcipinstance.h"

GLCIPInstance::GLCIPInstance(Digraph &g, DNodePosMap &posx, DNodePosMap &posy, ArcValueMap &influence, DNodeValueMap &threshold,
    DNodeListOfValuesMap &incentives, double alpha, int n, int m) : g(g), influence(influence), threshold(threshold),
    incentives(incentives), posx(posx), posy(posy)
{
    this->alpha = alpha;
    this->n = n;
    this->m = m;
}
