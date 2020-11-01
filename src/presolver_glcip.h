#ifndef __SCIP_PRESOLVER_GLCIP_H__
#define __SCIP_PRESOLVER_GLCIP_H__

#include "GLCIPBase.h"

class PresolverGLCIP : public ObjPresol
{
private:
    GLCIPInstance &instance;
    DNodeSCIPVarMap &x;
    ArcSCIPVarMap &z;
    DNodeInfSetsMap &infSet;

public:
    PresolverGLCIP(
        SCIP *scip,
        const char *p_name,        // name of branching rule
        GLCIPInstance &p_instance, // problem data
        DNodeSCIPVarMap &p_x,
        ArcSCIPVarMap &p_z,
        DNodeInfSetsMap &p_inf_set // influencing set data structure
    );
    virtual ~PresolverGLCIP();

    virtual SCIP_DECL_PRESOLEXEC(scip_exec);
};

#endif