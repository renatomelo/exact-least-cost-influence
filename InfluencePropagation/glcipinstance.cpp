#include "glcipinstance.h"

GLCIPInstance::GLCIPInstance(Digraph &g, DNodeStringMap &nodeName, DNodePosMap &posx, DNodePosMap &posy, ArcValueMap &influence,
    DNodeValueMap &threshold, DNodeValueVectorMap &incentives, double alpha, int n, int m) : g(g), influence(influence),
    threshold(threshold), incentives(incentives), nodeName(nodeName), posx(posx), posy(posy)
{
    this->alpha = alpha;
    this->n = n;
    this->m = m;
}
