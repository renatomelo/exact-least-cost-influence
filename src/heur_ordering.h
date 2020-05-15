#ifndef __HEURORDERING_H__
#define __HEURORDERING_H__

#include "GLCIPBase.h"

class HeurOrdering : public scip::ObjHeur
{
    GLCIPInstance &instance; /**< the instance for GLCIP */
    DNodeSCIPVarMap &x;
    ArcSCIPVarMap &z;
    DNodeSCIPVarsMap &xip;
    SCIP_SOL *sol; /**< current solution */

public:
    HeurOrdering(
        SCIP *scip,
        GLCIPBase &p_instance,
        DNodeSCIPVarMap &p_x,
        ArcSCIPVarMap &p_z,
        DNodeSCIPVarsMap &p_xip) : ObjHeur(scip,
                                      "MinInfluence",              //name
                                      "Ordering primal heuristic", //description
                                      'O',                         //display character of primal heuristic
                                      -1000000,                    //priority of the primal heuristic
                                      5,                           //frequency for calling primal heuristic
                                      0,                           //frequency offset for calling primal heuristic
                                      -1,                          //maximal depth level to call heuristic at (-1: no limit)
                                      SCIP_HEURTIMING_AFTERNODE,   //positions in the node solving loop where heuristic should be executed;
                                                                   //see definition of SCIP_HEURTIMING for possible values
                                      FALSE),                      //does the heuristic use a secondary SCIP instance?
                              instance(p_instance),
                              x(p_x),
                              z(p_z),
                              xip(p_xip),
                              sol(NULL)
    {
    }

    /** destructor */
    virtual ~HeurOrdering() {}

    /** destructor of primal heuristic to free user data (called when SCIP is exiting) */
    virtual SCIP_DECL_HEURFREE(scip_free);

    /** initialization method of primal heuristic (called after problem was transformed) */
    virtual SCIP_DECL_HEURINIT(scip_init);

    /** deinitialization method of primal heuristic (called before transformed problem is freed) */
    virtual SCIP_DECL_HEUREXIT(scip_exit);

    /** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  The primal heuristic may use this call to initialize its branch and bound specific data.
    *
    */
    virtual SCIP_DECL_HEURINITSOL(scip_initsol);

    /** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
    *
    *  This method is called before the branch and bound process is freed.
    *  The primal heuristic should use this call to clean up its branch and bound data.
    */
    virtual SCIP_DECL_HEUREXITSOL(scip_exitsol);

    /** execution method of primal heuristic
    *
    *  Searches for feasible primal solutions. The method is called in the node processing loop.
    *
    *  possible return values for *result:
    *  - SCIP_FOUNDSOL   : at least one feasible primal solution was found
    *  - SCIP_DIDNOTFIND : the heuristic searched, but did not find a feasible solution
    *  - SCIP_DIDNOTRUN  : the heuristic was skipped
    *  - SCIP_DELAYED    : the heuristic was skipped, but should be called again as soon as possible, disregarding
    *                      its frequency
    */
    virtual SCIP_DECL_HEUREXEC(scip_exec);
};

#endif