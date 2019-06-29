#ifndef __SCIP_EVENTHANDLER_GLCIP_H__
#define __SCIP_EVENTHANDLER_GLCIP_H__

#include "GLCIPBase.h"

class ObjEventhdlrGLCIP : public ObjEventhdlr
{
private:
    GLCIPInstance&          instance;

public:    
    //Constructs the event handler object for the branching rule
    ObjEventhdlrGLCIP(
        SCIP*               scip,
        const char*         p_name,
        GLCIPInstance&      p_instance
    );

    // Destructs the event handler
    virtual ~ObjEventhdlrGLCIP();

    virtual SCIP_DECL_EVENTFREE(scip_free);
    virtual SCIP_DECL_EVENTINIT(scip_init);
    virtual SCIP_DECL_EVENTEXIT(scip_exit);
    virtual SCIP_DECL_EVENTINITSOL(scip_initsol);
    virtual SCIP_DECL_EVENTEXITSOL(scip_exitsol);
    virtual SCIP_DECL_EVENTDELETE(scip_delete);
    virtual SCIP_DECL_EVENTEXEC(scip_exec);

    //SCIP_RETCODE SCIPincludeEventhdlrGLCIP(SCIP* scip);
};

#endif