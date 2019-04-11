#ifndef COLUMNGENERATOR_H
#define COLUMNGENERATOR_H

#include "mygraphlib.h"
#include <list>
#include <cassert>
#include "easyscip.h"
#include <utility>
#include <map>
#include <vector>
#include "star.h"
#include "glcipinstance.h"

class ColumnGenerator
{
private:
    static pair<double, vector<Arc>> solveKnapsack(vector<Arc> edges, vector<double> sizes, vector<double> costs, int capacity);

public:
    static Star getNewStar(GLCIPInstance instance, DNode root, ArcValueMap reducedCosts, double p);
};

#endif // KNAPSACKSOLVER_H
