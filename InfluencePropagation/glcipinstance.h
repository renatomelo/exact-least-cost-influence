#ifndef GLCIPINSTANCE_H
#define GLCIPINSTANCE_H

#include "mygraphlib.h"
#include <list>
#include <cassert>
#include "easyscip.h"
#include <utility>
#include <map>
#include <vector>
#include "star.h"

class GLCIPInstance
{
public:
    Digraph g;
    ArcValueMap influence;
    DNodeValueMap threshold;
    DNodeValueMap incentives;
    DNodeValueMap costs;
    double alpha;

    GLCIPInstance();
};

#endif // GLCIPINSTANCE_H
