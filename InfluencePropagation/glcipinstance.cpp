#include "glcipinstance.h"

GLCIPInstance::GLCIPInstance(Digraph g, ArcValueMap influence, DNodeValueMap threshold, DNodeValueMap incentives, DNodeValueMap costs,
    double alpha) : g(g), influence(influence), threshold(threshold), incentives(incentives), costs(costs)
{
    this->alpha = alpha;
}
