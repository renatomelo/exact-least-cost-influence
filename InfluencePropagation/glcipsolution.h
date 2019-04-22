#ifndef GLCIPSOLUTION_H
#define GLCIPSOLUTION_H

#include "mygraphlib.h"
#include <list>
#include <cassert>
#include "easyscip.h"
#include <utility>
#include <map>
#include <vector>
#include "star.h"

class GLCIPSolution
{
public:
    ArcBoolMap influence;
    DNodeValueMap incentives;

    GLCIPSolution(Digraph &g);
};

#endif // GLCIPSOLUTION_H
