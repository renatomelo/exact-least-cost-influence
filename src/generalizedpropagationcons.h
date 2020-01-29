/** 
 * Constraint handler for GLCIP generalized propagation constraints
 */
#ifndef __GENERALIZEDPROPAGATION_H__
#define __GENERALIZEDPROPAGATION_H__

#include "GLCIPBase.h"
#include "/opt/gurobi810/linux64/include/gurobi_c++.h"

class GeneralizedPropagation : public ObjConshdlr
{
private:
    GLCIPInstance &instance;
    DNodeSCIPVarMap &x;
    ArcSCIPVarMap &z;
    DNodeInfSetsMap &infSet;
    vector<Phi> &gpcrows;
    GRBEnv *env;
    //unsigned long int nGeneratedCons = 0;

public:
    GeneralizedPropagation(
        SCIP *scip,
        GLCIPInstance &p_instance,
        DNodeSCIPVarMap &p_x,
        ArcSCIPVarMap &p_z,
        DNodeInfSetsMap &p_infSet,
        vector<Phi> &p_gpcrows) : ObjConshdlr(scip,
                                              "GPC",
                                              "GLCIP generalized propagation constraints",
                                              1,     //priority for separation
                                              -2,    //priority for constraint enforcing
                                              1,     //priority for checking infeasibility (and propagation)
                                              1.0,   //frequency for separating cuts; zero means to separate only in the root node
                                              -1.0,  //frequency for propagating domains; zero means only preprocessing propagation
                                              1.0,   /* frequency for using all instead of only the useful constraints in separation,
                                                      *  propagation and enforcement, -1 for no eager evaluations, 0 for first only */
                                              0.0,   //maximal number of presolving rounds the constraint handler participates in (-1: no limit)
                                              FALSE, //should separation method be delayed, if other separators found cuts?
                                              FALSE,
                                              TRUE,
                                              SCIP_PROPTIMING_BEFORELP,
                                              SCIP_PRESOLTIMING_FAST),
                                  instance(p_instance),
                                  x(p_x),
                                  z(p_z),
                                  infSet(p_infSet),
                                  gpcrows(p_gpcrows)
    {
        env = new GRBEnv();
    }

    virtual ~GeneralizedPropagation() { delete env;}

    virtual SCIP_DECL_CONSDELETE(scip_delete);
    virtual SCIP_DECL_CONSTRANS(scip_trans);
    virtual SCIP_DECL_CONSSEPALP(scip_sepalp);
    virtual SCIP_DECL_CONSSEPASOL(scip_sepasol);
    virtual SCIP_DECL_CONSENFOLP(scip_enfolp);
    virtual SCIP_DECL_CONSENFOPS(scip_enfops);
    virtual SCIP_DECL_CONSCHECK(scip_check);
    virtual SCIP_DECL_CONSPROP(scip_prop);
    virtual SCIP_DECL_CONSLOCK(scip_lock);
    virtual SCIP_DECL_CONSDELVARS(scip_delvars);
    virtual SCIP_DECL_CONSPRINT(scip_print);

    /*  virtual SCIP_DECL_CONSHDLRISCLONEABLE(iscloneable)
    {
        return true;
    } */

    //virtual SCIP_DECL_CONSHDLRCLONE(ObjProbCloneable *clone);
    //virtual SCIP_DECL_CONSCOPY(scip_copy);

    SCIP_RETCODE greedSetExtensionHeur(
        SCIP *scip,
        SCIP_CONSHDLR *conshdlr, //the constraint handler itself
        SCIP_SOL *sol,           //primal solution that should be separated
        SCIP_RESULT *result      //pointer to store the result of the separation call
    );
    SCIP_RETCODE sepaGeneralizedPropCons(
        SCIP *scip,
        SCIP_CONSHDLR *conshdlr, //the constraint handler itself
        SCIP_SOL *sol,           //primal solution that should be separated
        SCIP_RESULT *result      //pointer to store the result of the separation call
                                 //        set<DNode> generalizedSet // set of vertices to be separated
    );

    SCIP_RETCODE exactSeparationGrbModel(
        SCIP *scip,
        SCIP_CONSHDLR *conshdlr, //the constraint handler itself
        SCIP_SOL *sol,           //primal solution that should be separated
        SCIP_RESULT *result     //pointer to store the result of the separation call
        );   

    SCIP_RETCODE exactSeparation(
        SCIP *scip,
        SCIP_CONSHDLR *conshdlr, //the constraint handler itself
        SCIP_SOL *sol,           //primal solution that should be separated
        SCIP_RESULT *result      //pointer to store the result of the separation call
    );

    SCIP_RETCODE createGenPropagationCons(
        SCIP *scip,
        SCIP_CONS **cons,
        const char *name);

    SCIP_RETCODE addGeneralizedPropCons(
        SCIP *scip,
        SCIP_CONSHDLR *conshdlr, //the constraint handler itself
        SCIP_SOL *sol,
        SCIP_RESULT *result,
        set<DNode> generalizedSet,
        DNode k,
        SCIP_Bool lifting);
};

#endif