#ifndef __SCIP_CONSHDLRARCMARKER_H__
#define __SCIP_CONSHDLRARCMARKER_H__

#include "GLCIPBase.h"
#include <string>

#define CONSHDLR_NAME "degree-constraint-handler"

namespace degreecons
{

    enum ConsType
    {
        UPPERBOUND = 1,
        LOWERBOUND = 0,
        FIXEDSIZE = 2
    };
    typedef enum ConsType CONSTYPE;

    class DegreeCons : public ObjConshdlr
    {
    private:
        GLCIPInstance &instance;
        ArcSCIPVarMap &z; /**< map of arc variables */
        DNodeSCIPVarMap &x;
        DNodeInfSetsMap &infSet;

    public:
        // default constructor
        DegreeCons(
            SCIP *scip,
            GLCIPInstance &p_instance,
            ArcSCIPVarMap &p_z, /**< map of arc variables */
            DNodeSCIPVarMap &p_x,
            DNodeInfSetsMap &p_infSet) : ObjConshdlr(scip,
                                                     CONSHDLR_NAME,                         //name
                                                     "stores the local branching decision", //description
                                                     0,                                     //separation priority
                                                     1,                                     //enforce priority
                                                     -1,                                     //checking feasibility priority
                                                     -1,                                    //separation frequency
                                                     10,                                     //propagatin frequeny
                                                     1,                                     //eager frequency
                                                     0,                                     //maximal presolving rounds
                                                     FALSE,                                 //delayed separation
                                                     FALSE,                                 //delayed propagation
                                                     TRUE,                                  //conshdlr should be skipped
                                                     SCIP_PROPTIMING_BEFORELP,              //propagation timing
                                                     SCIP_PRESOLTIMING_FAST),               //presolving timing
                                         instance(p_instance),
                                         z(p_z),
                                         x(p_x),
                                         infSet(p_infSet)
        {
        }

        // destructor
        virtual ~DegreeCons() {}

        virtual SCIP_DECL_CONSDELETE(scip_delete);
        virtual SCIP_DECL_CONSTRANS(scip_trans);
        //virtual SCIP_DECL_CONSSEPALP(scip_sepalp);
        //virtual SCIP_DECL_CONSSEPASOL(scip_sepasol);
        virtual SCIP_DECL_CONSENFOLP(scip_enfolp);
        virtual SCIP_DECL_CONSENFOPS(scip_enfops);
        virtual SCIP_DECL_CONSPROP(scip_prop);
        virtual SCIP_DECL_CONSACTIVE(scip_active);
        virtual SCIP_DECL_CONSDEACTIVE(scip_deactive);
        virtual SCIP_DECL_CONSLOCK(scip_lock);
        virtual SCIP_DECL_CONSDELVARS(scip_delvars);
        //virtual SCIP_DECL_CONSCOPY(scip_copy);
        virtual SCIP_DECL_CONSCHECK(scip_check);
        virtual SCIP_DECL_CONSPRINT(scip_print);

        CONSTYPE getTypeArcMarker(SCIP *scip, SCIP_CONS *cons);
        Arc getArc(SCIP *scip, SCIP_CONS *cons);
    };

    SCIP_RETCODE createDegreeCons(
        SCIP *scip,
        SCIP_CONS **cons,
        const char *name,
        DNode &vertex,
        int bound,
        CONSTYPE type,
        SCIP_NODE *node);

    SCIP_RETCODE createDegreeCons2(
        SCIP *scip,
        SCIP_CONS **cons,
        const char *name);

} // namespace degreecons
#endif