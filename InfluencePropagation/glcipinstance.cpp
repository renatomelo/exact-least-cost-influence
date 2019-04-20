#include "glcipinstance.h"

GLCIPInstance::GLCIPInstance(Digraph g, 
                            DNodePosMap posx,
                            DNodePosMap posy, 
                            ArcValueMap influence, 
                            DNodeValueMap threshold,
                            DNodeListOfValuesMap incentives, 
                            DNodeListOfValuesMap costs, 
                            double alpha, 
                            int n, 
                            int m) : g(g), 
                            influence(influence),
                            threshold(threshold), 
                            incentives(incentives), 
                            costs(costs), 
                            posx(posx), 
                            posy(posy)
{
    this->alpha = alpha;
}
