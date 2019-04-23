#ifndef GRAPHVIEWER_H
#define GRAPHVIEWER_H
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include "mygraphlib.h"
#include <string>
#include <vector>
#include <gurobi_c++.h>
#include "glcipinstance.h"
#include "glcipsolution.h"

class GraphViewer
{
    private:
        static vector<string> colors;

    public:
        static void ViewGLCIPSolution(GLCIPInstance &instance, GLCIPSolution &solution, string title);
};

#endif // GRAPHVIEWER_H
