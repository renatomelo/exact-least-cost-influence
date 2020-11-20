# InfluencePropagation
Branch-Cut-and-Price Algorithm for the Generalized Influence Propagation Problem

Dependences:
1. SCIP optimization framework version 6.*
2. Gurobi optimization solver version 8.1
3. Lemon graph library version 1.3

Compilation:

`cd exact-least-cost-influence/`

`make LPS=grb`

Execution:

`bin/glcip -i data/smallworld_10-8_p03.lcip -a arc -alpha 1`