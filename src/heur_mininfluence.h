#ifndef __HEURMININFLUENCE_H__
#define __HEURMININFLUENCE_H__

#include "GLCIPBase.h"

/** C++ wrapper */
class HeurMinInfluence : public scip::ObjHeur
{
   GLCIPInstance&     instance_;          /**< the underlying graph of the TSP */
   SCIP_SOL*          sol_;               /**< current solution */
   GRAPHEDGE**        tour_;              /**< tour induced by sol */

public:

   /** default constructor */
   HeurMinInfluence(
      SCIP* scip
      )
      : ObjHeur(scip, "2opt", "2-Opt heuristic for TSPs", 'K',-1000000, 1, 0, -1, SCIP_HEURTIMING_AFTERNODE, FALSE),
      graph_(0),
      ncalls_(0),
      sol_(NULL),
      tour_(NULL)
   {
   }

   /** destructor */
   virtual ~HeurMinInfluence()
   {
   } /*lint !e1540*/

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

   /** clone method which will be used to copy a objective plugin */
   virtual SCIP_DECL_HEURCLONE(ObjCloneable* clone); /*lint !e665*/

   /** returns whether the objective plugin is copyable */
   virtual SCIP_DECL_HEURISCLONEABLE(iscloneable)
   {
      return true;
   }
}; /*lint !e1712*/

#endif
