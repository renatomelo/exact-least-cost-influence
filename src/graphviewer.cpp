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

    GA.SetDefaultDNodeAttrib("color = Gray shape = circle style = bold fontsize = 18");

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
            //GA.SetColor(a, "Invis");
            GA.SetColor(a, "gray");
        }
    }

    GA.SetLabel(title);
    GA.View();
}

void GraphViewer::ViewGLCIPFracSolution(GLCIPInstance &instance, ArcValueMap &weight, string title){
    // assign names to the nodes
    DNodeStringMap nodeNames(instance.g);
    for(DNodeIt v(instance.g); v != INVALID; ++v){
        nodeNames[v] = "\"" + instance.nodeName[v] + "\"";
    }

    // set graph attributes for the visualizer
    DigraphAttributes GA(instance.g, nodeNames, instance.posx, instance.posy);

    GA.SetDefaultDNodeAttrib("color = Gray shape = ellipse style = bold fontsize = 20");

    // set arcs label according to weights
    for(ArcIt a(instance.g); a != INVALID; ++a){
        if(weight[a] > 0.00001){
            stringstream stream;
            stream << fixed << setprecision(4) << weight[a];
            string s = stream.str();

            GA.SetLabel(a, s);
        }
        else{
            GA.SetColor(a, "Invis");
        }
    }

    GA.SetLabel(title);
    GA.View();
}
