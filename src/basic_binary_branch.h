#ifndef __SCIP_BINARYBRANCH_GLCIP_H__
#define __SCIP_BINARYBRANCH_GLCIP_H__

#include "GLCIPBase.h"

class BasicBinaryBranch : public ObjBranchrule
{
private:
    GLCIPInstance&          instance;
    DNodeSCIPVarMap&        x;
    ArcSCIPVarMap&          z;
    //DNodeInfSetsMap&        infSet;

public:
    /** Constructs the branching rule object with the data needed */
    BasicBinaryBranch(
        SCIP*               scip,
        GLCIPInstance&      p_instance,     // problem data
        DNodeSCIPVarMap&    p_x,
        ArcSCIPVarMap&      p_z
        //DNodeInfSetsMap&    p_inf_set       // influencing set data structure 
    );

    // Destructs the branching rule object
    virtual ~BasicBinaryBranch();

    int leastFractional(
    SCIP *scip,
    SCIP_VAR **candidates,
    SCIP_Real *branchCandsFrac,
    int nCands);

    virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);
};

#endif