#include "GLCIPBase.h"

struct influencingSet{
    set<DNode> elements;
};

typedef Digraph::NodeMap< vector<influencingSet> > DNodeInfSetsMap;

class CovModelAllVariables: public GLCIPBase
{
public:
    CovModelAllVariables();
    ~CovModelAllVariables();
    static vector<influencingSet> powerSet(vector<DNode> neighbors);
    static double costInfluencingSet(GLCIPInstance &instance, DNode v, set<DNode> elements);
    static bool run(GLCIPInstance &instance, GLCIPSolution &solution, int timeLimit);
};
