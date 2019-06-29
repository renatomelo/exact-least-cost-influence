#include "event_glcip.h"
#include "pricer_glcip.h"

ObjEventhdlrGLCIP::ObjEventhdlrGLCIP(
    SCIP *scip,
    const char *p_name,
    GLCIPInstance &p_instance) : ObjEventhdlr(scip, p_name, "Event handler for update the pricer"),
                                 instance(p_instance)
{
}

ObjEventhdlrGLCIP::~ObjEventhdlrGLCIP() {}

/* SCIP_RETCODE ObjEventhdlrGLCIP::SCIPincludeEventhdlrGLCIP(SCIP* scip)
{
    return SCIP_OKAY;
} */

SCIP_DECL_EVENTFREE(ObjEventhdlrGLCIP::scip_free)
{
    return SCIP_OKAY;
}

SCIP_DECL_EVENTINIT(ObjEventhdlrGLCIP::scip_init)
{
    return SCIP_OKAY;
}

SCIP_DECL_EVENTEXIT(ObjEventhdlrGLCIP::scip_exit)
{
    return SCIP_OKAY;
}

SCIP_DECL_EVENTINITSOL(ObjEventhdlrGLCIP::scip_initsol)
{
    SCIP_VAR **vars = SCIPgetVars(scip);
    int nvars = SCIPgetNVars(scip);

    for (int i = 0; i < nvars; i++)
    {
        SCIP_CALL(SCIPcatchVarEvent(scip, 
                                    vars[i],
                                    SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_BOUNDCHANGED, 
                                    eventhdlr, 
                                    NULL, NULL));
    }

    return SCIP_OKAY;
}

SCIP_DECL_EVENTEXITSOL(ObjEventhdlrGLCIP::scip_exitsol)
{
    //SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, -1));
    SCIP_VAR **vars = SCIPgetVars(scip);
    int nvars = SCIPgetNVars(scip);

    for (int i = 0; i < nvars; i++)
    {
        SCIP_CALL( SCIPdropVarEvent(scip, 
                                    vars[i],
                                    SCIP_EVENTTYPE_GBDCHANGED | SCIP_EVENTTYPE_BOUNDCHANGED, 
                                    eventhdlr, 
                                    NULL, 
                                    -1) );
    }
    return SCIP_OKAY;
}

SCIP_DECL_EVENTDELETE(ObjEventhdlrGLCIP::scip_delete)
{
    return SCIP_OKAY;
}

SCIP_DECL_EVENTEXEC(ObjEventhdlrGLCIP::scip_exec)
{
    /* SCIP_SOL* sol = SCIPgetBestSol(scip);    
    assert(sol != NULL);
    assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND);
    SCIP_Real val = SCIPgetSolOrigObj(scip, sol);
    SCIPinfoMessage(scip, NULL, "found new best solution with solution value <%g>\n", val); */

    SCIP_VAR *var;
    ObjPricerGLCIP *pricer;
    SCIP_PRICERDATA *pricerData;
    SCIP_VAR **origninalVars;

    pricer = static_cast<ObjPricerGLCIP *>(SCIPfindObjPricer(scip, "GLCIP_pricer"));
    assert(pricer != NULL);

    //TODO maybe define the pricer data structure
    //pricerData = pricer->getPricerData();

    assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARDELETED);
    var = SCIPeventGetVar(event);
    assert(var != NULL);

    /* SCIP_EVENTTYPE_VARFIXED
    SCIP_EVENTTYPE_OBJCHANGED
    SCIP_EVENTTYPE_NODEFOCUSED */
    
    // adaptation
    origninalVars = SCIPgetOrigVars(scip);
    assert(origninalVars != NULL);

    for (int i = 0; i < SCIPgetNOrigVars(scip); i++)
    {
        SCIPinfoMessage(scip, NULL, "deleted var <%s>\n", SCIPvarGetName(origninalVars[i]));
    }
    
    return SCIP_OKAY;
}