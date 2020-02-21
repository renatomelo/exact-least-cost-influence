#ifndef __HEURDUALBOUND_H__
#define __HEURDUALBOUND_H__

#include "GLCIPBase.h"

class HeurDualBound : public scip::ObjRelax
{
    GLCIPInstance &instance; /**< the instance for GLCIP */
    DNodeSCIPVarMap &x;
    ArcSCIPVarMap &z;
    DNodeSCIPVarsMap &xip;
    SCIP_SOL *sol_; /**< current solution */

public:
    // default constructor
    HeurDualBound(
        SCIP *scip,
        GLCIPInstance &p_instance,
        DNodeSCIPVarMap &p_x,
        ArcSCIPVarMap &p_z,
        DNodeSCIPVarsMap &p_xip) : ObjRelax(scip,
                                       "heuristic-dual-bound",
                                       "Heuristic dual bound for GLCIP",
                                       -1.0, //priority of the relaxator (negative: after LP, non-negative: before LP)
                                       5,    //frequency for calling relaxator
                                       0),   //Does the relaxator contain all cuts in the LP?
                              instance(p_instance),
                              x(p_x),
                              z(p_z),
                              xip(p_xip),
                              sol_(NULL)
    {
    }
    //destructor
    virtual ~HeurDualBound() {}

    virtual SCIP_DECL_RELAXFREE(scip_free);
    virtual SCIP_DECL_RELAXINIT(scip_init);
    virtual SCIP_DECL_RELAXEXIT(scip_exit);
    virtual SCIP_DECL_RELAXINITSOL(scip_initsol);
    virtual SCIP_DECL_RELAXEXITSOL(scip_exitsol);

    /** execution method of relaxator
     *
     *  The method is called in the node processing loop. It solves the current subproblem's relaxation.
     *  Like the LP relaxation, the relaxator should only operate on COLUMN variables.
     *
     *  input:
     *  - scip            : SCIP main data structure
     *  - relax           : the relaxator itself
     *  - lowerbound      : pointer to store a lowerbound for the current node
     *  - result          : pointer to store the result of the relaxation call
     *
     *  possible return values for *result (if more than one applies, the first in the list should be used):
     *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
     *  - SCIP_CONSADDED  : an additional constraint was generated, and the relaxator should not be called again on the
     *                      same relaxation
     *  - SCIP_REDUCEDDOM : a variable's domain was reduced, and the relaxator should not be called again on the same
     *                      relaxation
     *  - SCIP_SEPARATED  : a cutting plane was generated, and the relaxator should not be called again on the same relaxation
     *  - SCIP_SUCCESS    : the relaxator solved the relaxation and should not be called again on the same relaxation
     *  - SCIP_SUSPENDED  : the relaxator interrupted its solving process to wait for additional input (e.g. cutting
     *                      planes); however, it is able to continue the solving in order to improve the dual bound
     *  - SCIP_DIDNOTRUN  : the relaxator was skipped
     */
    virtual SCIP_DECL_RELAXEXEC(scip_exec);
};

#endif