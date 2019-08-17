#ifndef __SCIP_BRANCH_GLCIP_H__
#define __SCIP_BRANCH_GLCIP_H__

#include "GLCIPBase.h"
#include "consarcmarker.h"

class ObjBranchruleGLCIP : public ObjBranchrule
{
private:
    GLCIPInstance&          instance;
    DNodeSCIPVarMap&        x;
    ArcSCIPVarMap&          z;
    DNodeInfSetsMap&        infSet;
    ArcBoolMap&             isAble;

public:
    /** Constructs the branching rule object with the data needed */
    ObjBranchruleGLCIP(
        SCIP*               scip,
        const char*         p_name,         // name of branching rule
        GLCIPInstance&      p_instance,     // problem data
        DNodeSCIPVarMap&    p_x,
        ArcSCIPVarMap&      p_z,
        DNodeInfSetsMap&    p_inf_set,       // influencing set data structure 
        ArcBoolMap&         p_isAble
    );

    // Destructs the branching rule object
    virtual ~ObjBranchruleGLCIP();

    //virtual SCIP_DECL_BRANCHINIT(scip_init);
    virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);
};

#endif