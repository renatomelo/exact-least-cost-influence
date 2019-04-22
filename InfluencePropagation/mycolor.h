// =============================================================
//
//  Routines for color codification
//
//  Send comments/corrections to Flavio K. Miyazawa.
//
// =============================================================
#ifndef MYCOLOR_DEFINE
#define MYCOLOR_DEFINE
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include "myutils.h"

using namespace std;

typedef struct {
  string colorname, rgb;
  int rgb_red,rgb_green,rgb_blue;
  int nextvisualdistinctcolor;
} colorelement;

class ColorClass {
public:
  ColorClass(void); // Use colors defined in the routine init
  // ColorClass(string filename);//Add color from a file (also includes the colors in init)
  // some functions
  string name(int colorcode);
  int insert(string colorname, string rgb);
  string rgb(int colorcode);
  bool exists(string colorname);
  int size();
  void init();
  void print();
  void list();
  // dictionary that returns the code of a given color name
  std::map<string,int> code;
  vector<colorelement> vecolor;
  vector<int> ith_visualdistinctcolor;
  // int MaxColors;
  int firstvisualdistinctcolor,lastvisualdistinctcolor,
    quantvisualdistinctcolor;
};
extern ColorClass ColorG;

inline int ColorCode(string colorname){return(ColorG.code[colorname]);}
inline int ColorSize(){return(ColorG.vecolor.size());}
inline string ColorRGB(string colorname){return(ColorG.rgb(ColorG.code[colorname]));}
inline string ColorRGB(int  colorcode){return(ColorG.rgb(colorcode));}
inline bool ColorExists(string colorname){return(ColorG.exists(colorname));}
inline string ColorName(int code){return(ColorG.name(code));}
inline int ColorInsert(string colorname,string rgb){return(ColorG.insert(colorname,rgb));}
inline void ColorPrint(){ColorG.print();}
inline void ColorList(){ColorG.list();}
int InsertVisualDistinctColor(string colorname);
int NextVisualDistinctColor(int colorcode);
int FirstVisualDistinctColor();
int SizeVisualDistinctColor();
int ith_VisualDistinctColor(int i);
inline string ith_VisualDistinctColorName(int i){return(ColorG.name(ith_VisualDistinctColor(i)));}
#endif
