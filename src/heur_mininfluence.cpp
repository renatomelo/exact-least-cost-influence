#include "heur_mininfluence.h"

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_DECL_HEURFREE(HeurMinInfluence::scip_free)
{
   return SCIP_OKAY;
} 


/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_DECL_HEURINIT(HeurMinInfluence::scip_init)
{
   cout << "SCIP_DECL_HEURINIT\n";
   return SCIP_OKAY;
} 


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_DECL_HEUREXIT(HeurMinInfluence::scip_exit)
{
   return SCIP_OKAY;
} 


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 */
SCIP_DECL_HEURINITSOL(HeurMinInfluence::scip_initsol)
{
   cout << "SCIP_DECL_HEURINITSOL\n";
   return SCIP_OKAY;
} 

   
/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_DECL_HEUREXITSOL(HeurMinInfluence::scip_exitsol)
{
   cout << "SCIP_DECL_HEUREXITSOL\n";
   return SCIP_OKAY;
} 


/** execution method of primal heuristic 2-Opt */
SCIP_DECL_HEUREXEC(HeurMinInfluence::scip_exec)
{  
   cout << "SCIP_DECL_HEUREXEC\n";
   assert( heur != NULL );
   SCIP_SOL* sol = SCIPgetBestSol( scip );
   bool newsol;

   // check whether a new solution was found meanwhile
   if( sol != sol_ )
   {
      sol_ = sol;
      newsol = true;
   }
   else
      newsol = false;

   // get tour from sol and sort edges by length, if new solution was found
   if( newsol )
   {
      //todo

   }


   return SCIP_OKAY;
} 

/** clone method which will be used to copy a objective plugin */
/* SCIP_DECL_HEURCLONE(scip::ObjCloneable* HeurMinInfluence::clone)
{
   return new HeurMinInfluence(scip);
} */
