#include "graphviewer.h"

// just some colors that we can use
vector<string> GraphViewer::colors = {"Gray", "Red", "Blue", "Magenta", "Orange", "Green", "Pink", "Brown", "Purple", "Yellow",
                            "Cyan", "Olive", "Beige", "Lime", "Mint", "Teal", "Navy", "Lavender"};

void GraphViewer::ViewGLCIPSolution(GLCIPInstance &instance, GLCIPSolution &solution, string title){
    // assign names to the nodes
    DNodeStringMap nodeNames(instance.g);
    for(DNodeIt v(instance.g); v != INVALID; ++v){
        stringstream stream;
        stream << fixed << setprecision(2) << solution.incentives[v];
        string s = stream.str();

        nodeNames[v] = "\"" + instance.nodeName[v] + "(" + s + ")" + "\"";
    }

    // set graph attributes for the visualizer
    DigraphAttributes GA(instance.g, nodeNames, instance.posx, instance.posy);

    GA.SetDefaultDNodeAttrib("color = Gray shape = ellipse style = bold fontsize = 20");

    // color nodes according to the solution
    for(DNodeIt v(instance.g); v != INVALID; ++v){
        if(solution.incentives[v] > 0)
            GA.SetColor(v, GraphViewer::colors[1]);
        else
            GA.SetColor(v, GraphViewer::colors[0]);
    }

    // remove arcs that are not in the solution
    for(ArcIt a(instance.g); a != INVALID; ++a){
        if(!solution.influence[a]){
            GA.SetColor(a, "Invis");
        }
    }

    GA.SetLabel(title);
    GA.View();
}
