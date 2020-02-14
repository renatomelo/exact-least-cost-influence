#include "heur_dualbound.h"

SCIP_DECL_RELAXFREE(HeurDualBound::scip_free)
{
    //cout << "RELAXFREE()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXINIT(HeurDualBound::scip_init)
{
    //cout << "RELAXEINIT()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXEXIT(HeurDualBound::scip_exit)
{
    //cout << "RELAXEXIT()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXINITSOL(HeurDualBound::scip_initsol)
{
    //cout << "RELAXINITSOL()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXEXITSOL(HeurDualBound::scip_exitsol)
{
    //cout << "RELAXEXITSOL()" << endl;
    return SCIP_OKAY;
}

SCIP_DECL_RELAXEXEC(HeurDualBound::scip_exec)
{
    cout << "RELAXEXEC()" << endl;

    SCIP_CONS **conss;
    int nconss;
    SCIP_Real relaxval;
    SCIP_Bool valid;

    conss = SCIPgetConss(scip);
    nconss = SCIPgetNConss(scip);

    for (int c = 0; c < nconss; ++c)
    {
        const char *conshdlrname;

        conshdlrname = SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c]));

        cout << conshdlrname << endl;
    }

    *result = SCIP_DIDNOTRUN;
    *lowerbound = -SCIPinfinity(scip);

    //call the heuristic here
    relaxval = getLowerBound(scip, instance, x, z);

    printf("relaxation bound = %e \n", relaxval);

    //if( SCIPgetStatus(relaxscip) == SCIP_STATUS_OPTIMAL )
    {
        /* store relaxation solution in original SCIP if it improves the best relaxation solution thus far */
        if ((!SCIPisRelaxSolValid(scip)) || SCIPisGT(scip, relaxval, SCIPgetRelaxSolObj(scip)))
        {
            cout << "Setting LP relaxation solution, which improved upon earlier solution\n";
            SCIP_CALL(SCIPclearRelaxSolVals(scip));

            //it is necessary to set the new solution? I think so
            for (i = 0; i < SCIPgetNVars(scip); ++i)
            {
                SCIP_VAR *relaxvar;
                SCIP_Real solval;

                relaxvar = SCIPhashmapGetImage(varmap, SCIPgetVars(scip)[i]);
                assert(relaxvar != NULL);

                solval = SCIPgetSolVal(relaxscip, SCIPgetBestSol(relaxscip), relaxvar);

                SCIP_CALL(SCIPsetRelaxSolVal(scip, SCIPgetVars(scip)[i], solval));
            }

            //TODO: if we found a strongly connected component add a cuting plane

            /* mark relaxation solution to be valid and inform SCIP that the relaxation included all LP rows */
            SCIP_CALL(SCIPmarkRelaxSolValid(scip, TRUE));
        }

        SCIPdebugMsg(scip, "LP lower bound = %g\n", relaxval);
        *lowerbound = relaxval;
        *result = SCIP_SUCCESS;
    }

    return SCIP_OKAY;
}