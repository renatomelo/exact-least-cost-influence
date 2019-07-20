#ifndef __SCIP_CONSHDLRARCMARKER_H__
#define __SCIP_CONSHDLRARCMARKER_H__

#include "GLCIPBase.h"

#define CONSHDLR_NAME "arcmarker"

namespace arcmarker
{

enum ConsType
{
    WITH = 1,
    WITHOUT = 0
};
typedef enum ConsType CONSTYPE;

class ConshdlrArcMarker : public ObjConshdlr
{
private:
    GLCIPInstance &instance;
    ArcSCIPVarMap &z; /**< map of arc variables */

public:
    // default constructor
    ConshdlrArcMarker(
        SCIP *scip,
        GLCIPInstance &p_instance,
        ArcSCIPVarMap &p_z /**< map of arc variables */
        ) : ObjConshdlr(scip,
                        CONSHDLR_NAME,                         //name
                        "stores the local branching decision", //description
                        0,                                     //separation priority
                        0,                                     //enforce priority
                        9999999,                               //checking feasibility priority
                        0,                                     //separation frequency
                        1,                                     //propagatin frequeny
                        1,                                     //eager frequency
                        0,                                     //maximal presolving rounds
                        FALSE,                                 //delayed separation
                        FALSE,                                 //delayed propagation
                        TRUE,                                  //conshdlr should be skipped
                        SCIP_PROPTIMING_BEFORELP,              //propagation timing
                        SCIP_PRESOLTIMING_FAST),               //presolving timing
            instance(p_instance),
            z(p_z)
    {
    }

    // destructor
    virtual ~ConshdlrArcMarker() {}

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

    CONSTYPE getTypeArcMarker(SCIP *scip, SCIP_CONS *cons);
    Arc getArc(SCIP *scip, SCIP_CONS *cons);
};

SCIP_RETCODE createConsArcMarker(
    SCIP *scip,
    SCIP_CONS **cons,
    const char *name,
    SCIP_VAR* arcVar,
    ArcIt& arc,
    CONSTYPE type,
    SCIP_NODE *node);

} // namespace arcmarker
#endif