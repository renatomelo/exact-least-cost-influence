#ifndef __SCIP_BRANCH_GLCIP_H__
#define __SCIP_BRANCH_GLCIP_H__

#include "GLCIPBase.h"

class ObjBranchruleGLCIP : public ObjBranchrule
{
private:
    GLCIPInstance&          instance;
    DNodeInfSetsMap&        infSet;

public:
    /** Constructs the branching rule object with the data needed */
    ObjBranchruleGLCIP(
        SCIP*               scip,
        const char*         p_name,         // name of branching rule
        GLCIPInstance&      p_instance,     // problem data
        DNodeInfSetsMap&    p_inf_set       // influencing set data structure 
    );

    // Destructs the branching rule object
    virtual ~ObjBranchruleGLCIP();

    //virtual SCIP_DECL_BRANCHINIT(scip_init);
    virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);

};

#endif