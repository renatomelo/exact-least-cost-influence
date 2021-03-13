# Branch-and-Cut-and-Price for Generalized Least Cost Influence Problem
For study purposes, this project implements all the mathematical models and algorithms proposed by Fischetti et al. (2018) in the following paper.

- [Least cost influence propagation in (social) networks](https://link.springer.com/article/10.1007/s10107-018-1288-y)

Mathematical models:
1. Coverage Model (branch-and-cut-and-price)
2. Arc model (branch-and-cut)

Callbacks:
1. Cycle elimination constraints (separation problem)
    - Exact separation
2. Generalized propagation constraints (separation problem)
    - Exact and heuristic separation
3. Pricing problem for column generation (dynamic programming algorithm)

In addition to the algorithms and models of the paper above, this project contains additional ideas of algorithms to solve the problem, such as:
1. Branching rules
2. Presolving methods
3. Dual bond algorithm
4. Primal heuristics

## Dependences:
1. SCIP optimization framework version 6.*
2. Gurobi optimization solver version 8.1
3. Lemon graph library version 1.3

## Compilation:

`cd exact-least-cost-influence/`

`make LPS=grb`

Execution:

`bin/glcip -i data/smallworld_10-8_p03.lcip -a covcg -alpha 1`
