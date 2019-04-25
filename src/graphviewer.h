#ifndef GRAPHVIEWER_H
#define GRAPHVIEWER_H
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include "mygraphlib.h"
#include <string>
#include <vector>
#include "glcipinstance.h"
#include "glcipsolution.h"

class GraphViewer
{
    private:
        static vector<string> colors;

    public:
        static void ViewGLCIPSolution(GLCIPInstance &instance, GLCIPSolution &solution, string title);
        static void ViewGLCIPFracSolution(GLCIPInstance &instance, ArcValueMap &weight, string title);
};

#endif // GRAPHVIEWER_H
