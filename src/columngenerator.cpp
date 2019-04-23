#include "columngenerator.h"

/* Recall that the pricing subproblem is a knapsack problem in its minimization form.
 * The algorithm is just as the traditional knapsack algorithm: we use a dynamic programming
 * matrix DP[i][j], where DP[i][j] indicates an optimal solution that fills at least j with the items
 * 1, ..., i. From this, we can build the algorithm by induction.
 *
 * Base Case:
 * DP[0][j] = 0, if j == 0
 *          = infty, otherwise
 *
 * Induction Step: let s_i and c_i be the size and cost of item i, respectively. Then:
 * DP[i][j] = min(DP[i - 1][j - s_i] + c_i, DP[i - 1][j])
 * if j - s_i < 0, then we just use DP[i - 1][j]
 *
 * Since DP[i][j] just uses DP[i - 1], we can do this with only two vectors, instead of using a matrix
 */
pair<double, vector<Arc>> ColumnGenerator::solveKnapsack(vector<Arc> edges, vector<double> sizes, vector<double> costs, int capacity){
    std::vector<pair<double, vector<Arc>>> DP1 = new std::vector<pair<double, vector<Arc>>>(capacity);
    std::vector<pair<double, vector<Arc>>> DP2 = new std::vector<pair<double, vector<Arc>>>(capacity);

    // base case
    DP1[0] = 0;
    for(int j = 1; j <= capacity; j++){
        DP1[j].first = 10000;
    }

    // induction step
    for(int i = 1; i <= edges.size(); i++){
        for(int j = 0; j <= capacity; j++){
            double a = 10000;
            if(j - sizes[i] >= 0){
                a = DP1[j - sizes[i]] + costs[i];
            }
            double b = DP1[j];

            if(a < b){
                DP2[j] = DP1[j - sizes[i]];
                DP2[j].first += costs[i];
                DP2[j].second.push_back(edges[i]);
            }
            else{
                DP2[j] = DP1[j];
            }
        }
        DP1 = DP2;
    }

    // return the answer
    return DP2[capacity];
}

Star ColumnGenerator::getNewStar(GLCIPInstance instance, DNode root, ArcValueMap reducedCosts, double p){
    double aux_cap = p;
    double aux_cost = 0;
    vector<Arc> edges;
    vector<double> sizes;
    vector<double> costs;
    Star star;
    star.weight = p;
    star.root = root;

    // put arcs with negative reduced costs in the star
    for(OutArcIt a(instance.g, root); a != INVALID; ++a){
        if(reducedCosts[a] < 0){
            star.edges.push_back(a);
            aux_cap += instance.influence[a];
            aux_cost += reducedCosts[a];
        }
        else{
            edges.push_back(a);
            sizes.push_back(instance.influence[a]);
            costs.push_back(reducedCosts[a]);
        }
    }

    // solve the knapsack
    pair<double, vector<Arc>> sol = ColumnGenerator::solveKnapsack(edges, sizes, costs, capacity - aux_cap);

    // if best solution has negative cost, return this star to be added to the model
    if(sol.first + aux_cost < 0){
        for(auto a : sol.second){
            star.edges.push_back(a);
        }

        return star;
    }

    // otherwise, no star rooted at root can currently improve the model
    else{
        return NULL;
    }
}
