#ifndef ARCMODEL_H
#define ARCMODEL_H

#include "GLCIPBase.h"

class ArcModel: public GLCIPBase
{
    public:
    //       static void addCuttingPlanes(SCIP *scip, GLCIPInstance &instance, DNodeSCIPVarMap &x, ArcSCIPVarMap &z);
        static bool run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit);
};

#endif // ARCMODEL_H
