#ifndef STAR_H
#define STAR_H

#include "mygraphlib.h"
#include <list>
#include <cassert>
#include <utility>
#include <set>
#include <stack>
#include <vector>

class Star{
public:
    vector<Arc> edges;
    DNode root;
    double weight;
};

#endif // STAR_H
