//
// Class to use color names and codes.
//
// In the future, use this structure to represent colors in mygraphlib.
//
// Send comments/corrections to Flavio K. Miyazawa
//
#include "myutils.h"
#include "mycolor.h"

ColorClass ColorG; // As a global variable

int color2int(colorelement c){
  return(c.rgb_red*255*255+c.rgb_green*255+c.rgb_blue);}
bool compare_color(const colorelement &lhs, const colorelement &rhs){
  return color2int(lhs) < color2int(rhs);}

string ColorClass::name(int colorcode){return(this->vecolor[colorcode].colorname);}
string ColorClass::rgb(int colorcode){return(this->vecolor[colorcode].rgb);}
int ColorClass::size(){return(this->vecolor.size());}
bool ColorClass::exists(string colorname){
  if(this->code.find(colorname)!=this->code.end())
    return(true);else return(false);}
int ColorClass::insert(string colorname, string rgb){
  int pos;
  if (this->exists(colorname)) {
    pos=this->code[colorname];
    this->vecolor[pos].rgb = rgb;
    this->vecolor[pos].rgb_red = hex2int(rgb.substr(1,2));
    this->vecolor[pos].rgb_green = hex2int(rgb.substr(3,2));
    this->vecolor[pos].rgb_blue = hex2int(rgb.substr(5,2));
  }else{
    colorelement e;
    e.colorname = colorname;
    e.rgb = rgb;
    e.rgb_red = hex2int(rgb.substr(1,2));
    e.rgb_green = hex2int(rgb.substr(3,2));
    e.rgb_blue = hex2int(rgb.substr(5,2));
    e.nextvisualdistinctcolor = 0;
    this->vecolor.push_back(e);
    pos = this->vecolor.size()-1;// position of the last element
    this->code[colorname] = pos; }
  return(pos);
}

// you must not insert a same color more than once
int InsertVisualDistinctColor(string colorname){
  int colorcode = ColorCode(colorname);
  if (ColorG.firstvisualdistinctcolor==-1){
    ColorG.firstvisualdistinctcolor = colorcode;
    ColorG.lastvisualdistinctcolor = colorcode;
    ColorG.vecolor[colorcode].nextvisualdistinctcolor = colorcode;
  }else{
    ColorG.vecolor[colorcode].nextvisualdistinctcolor = ColorG.firstvisualdistinctcolor;
    ColorG.vecolor[ColorG.lastvisualdistinctcolor].nextvisualdistinctcolor = colorcode;
    ColorG.lastvisualdistinctcolor = colorcode;
  }
  ColorG.ith_visualdistinctcolor.push_back(colorcode);
  ColorG.quantvisualdistinctcolor++;
  return(colorcode);
}
int NextVisualDistinctColor(int colorcode){
  if (ColorG.firstvisualdistinctcolor==-1) return(-1);
  return(ColorG.vecolor[colorcode].nextvisualdistinctcolor);
}
int FirstVisualDistinctColor(){return(ColorG.firstvisualdistinctcolor);}
int SizeVisualDistinctColor(){return(ColorG.quantvisualdistinctcolor);}

int ith_VisualDistinctColor(int i){
  return(ColorG.ith_visualdistinctcolor[i]);}

void ColorClass::init(void)
{
  this->firstvisualdistinctcolor = -1;
  this->lastvisualdistinctcolor = -1;
  this->quantvisualdistinctcolor = 0;
  this->code["NoColor"] = -1; // This is a special code to represent no given color
  this->insert("Aqua" , "#00FFFF");// In some cases, the color 0 (first color) is default
  this->insert("Cyan" , "#00FFFF");
  this->insert("Gray" , "#808080");
  this->insert("Grey" , "#808080");
  this->insert("Black" , "#000000");
  this->insert("Navy" , "#000080");
  this->insert("DarkBlue" , "#00008B");
  this->insert("MediumBlue" , "#0000CD");
  this->insert("Blue" , "#0000FF");
  this->insert("DarkGreen" , "#006400");
  this->insert("Green" , "#008000");
  this->insert("Teal" , "#008080");
  this->insert("DarkCyan" , "#008B8B");
  this->insert("DeepSkyBlue" , "#00BFFF");
  this->insert("DarkTurquoise" , "#00CED1");
  this->insert("MediumSpringGreen" , "#00FA9A");
  this->insert("Lime" , "#00FF00");
  this->insert("SpringGreen" , "#00FF7F");
  this->insert("MidnightBlue" , "#191970");
  this->insert("DodgerBlue" , "#1E90FF");
  this->insert("LightSeaGreen" , "#20B2AA");
  this->insert("ForestGreen" , "#228B22");
  this->insert("SeaGreen" , "#2E8B57");
  this->insert("DarkSlateGray" , "#2F4F4F");
  this->insert("DarkSlateGrey" , "#2F4F4F");
  this->insert("LimeGreen" , "#32CD32");
  this->insert("MediumSeaGreen" , "#3CB371");
  this->insert("Turquoise" , "#40E0D0");
  this->insert("RoyalBlue" , "#4169E1");
  this->insert("SteelBlue" , "#4682B4");
  this->insert("DarkSlateBlue" , "#483D8B");
  this->insert("MediumTurquoise" , "#48D1CC");
  this->insert("Indigo" , "#4B0082");
  this->insert("DarkOliveGreen" , "#556B2F");
  this->insert("CadetBlue" , "#5F9EA0");
  this->insert("CornflowerBlue" , "#6495ED");
  this->insert("RebeccaPurple" , "#663399");
  this->insert("MediumAquaMarine" , "#66CDAA");
  this->insert("DimGrey" , "#696969");
  this->insert("DimGray" , "#696969");
  this->insert("SlateBlue" , "#6A5ACD");
  this->insert("OliveDrab" , "#6B8E23");
  this->insert("SlateGrey" , "#708090");
  this->insert("SlateGray" , "#708090");
  this->insert("LightSlateGray" , "#778899");
  this->insert("LightSlateGrey" , "#778899");
  this->insert("MediumSlateBlue" , "#7B68EE");
  this->insert("LawnGreen" , "#7CFC00");
  this->insert("Maroon" , "#800000");
  this->insert("Chartreuse" , "#7FFF00");
  this->insert("Purple" , "#800080");
  this->insert("Aquamarine" , "#7FFFD4");
  this->insert("Olive" , "#808000");
  this->insert("SkyBlue" , "#87CEEB");
  this->insert("LightSkyBlue" , "#87CEFA");
  this->insert("BlueViolet" , "#8A2BE2");
  this->insert("DarkRed" , "#8B0000");
  this->insert("DarkMagenta" , "#8B008B");
  this->insert("SaddleBrown" , "#8B4513");
  this->insert("DarkSeaGreen" , "#8FBC8F");
  this->insert("LightGreen" , "#90EE90");
  this->insert("MediumPurple" , "#9370DB");
  this->insert("DarkViolet" , "#9400D3");
  this->insert("PaleGreen" , "#98FB98");
  this->insert("DarkOrchid" , "#9932CC");
  this->insert("YellowGreen" , "#9ACD32");
  this->insert("Sienna" , "#A0522D");
  this->insert("Brown" , "#A52A2A");
  this->insert("DarkGray" , "#A9A9A9");
  this->insert("DarkGrey" , "#A9A9A9");
  this->insert("LightBlue" , "#ADD8E6");
  this->insert("GreenYellow" , "#ADFF2F");
  this->insert("PaleTurquoise" , "#AFEEEE");
  this->insert("LightSteelBlue" , "#B0C4DE");
  this->insert("PowderBlue" , "#B0E0E6");
  this->insert("FireBrick" , "#B22222");
  this->insert("DarkGoldenRod" , "#B8860B");
  this->insert("MediumOrchid" , "#BA55D3");
  this->insert("RosyBrown" , "#BC8F8F");
  this->insert("DarkKhaki" , "#BDB76B");
  this->insert("Silver" , "#C0C0C0");
  this->insert("LightGray" , "#D3D3D3");
  this->insert("LightGrey" , "#D3D3D3");
  this->insert("MediumVioletRed" , "#C71585");
  this->insert("IndianRed" , "#CD5C5C");
  this->insert("Peru" , "#CD853F");
  this->insert("Chocolate" , "#D2691E");
  this->insert("Tan" , "#D2B48C");
  this->insert("Thistle" , "#D8BFD8");
  this->insert("Orchid" , "#DA70D6");
  this->insert("GoldenRod" , "#DAA520");
  this->insert("PaleVioletRed" , "#DB7093");
  this->insert("Crimson" , "#DC143C");
  this->insert("Gainsboro" , "#DCDCDC");
  this->insert("Plum" , "#DDA0DD");
  this->insert("BurlyWood" , "#DEB887");
  this->insert("LightCyan" , "#E0FFFF");
  this->insert("Lavender" , "#E6E6FA");
  this->insert("DarkSalmon" , "#E9967A");
  this->insert("Violet" , "#EE82EE");
  this->insert("PaleGoldenRod" , "#EEE8AA");
  this->insert("LightCoral" , "#F08080");
  this->insert("Khaki" , "#F0E68C");
  this->insert("AliceBlue" , "#F0F8FF");
  this->insert("HoneyDew" , "#F0FFF0");
  this->insert("Azure" , "#F0FFFF");
  this->insert("SandyBrown" , "#F4A460");
  this->insert("Wheat" , "#F5DEB3");
  this->insert("Beige" , "#F5F5DC");
  this->insert("WhiteSmoke" , "#F5F5F5");
  this->insert("MintCream" , "#F5FFFA");
  this->insert("Mint" , "#AAFFC3");
  this->insert("GhostWhite" , "#F8F8FF");
  this->insert("Salmon" , "#FA8072");
  this->insert("AntiqueWhite" , "#FAEBD7");
  this->insert("Linen" , "#FAF0E6");
  this->insert("LightGoldenRodYellow" , "#FAFAD2");
  this->insert("OldLace" , "#FDF5E6");
  this->insert("Red" , "#FF0000");
  this->insert("Fuchsia" , "#FF00FF");
  this->insert("Magenta" , "#FF00FF");
  this->insert("DeepPink" , "#FF1493");
  this->insert("OrangeRed" , "#FF4500");
  this->insert("Tomato" , "#FF6347");
  this->insert("HotPink" , "#FF69B4");
  this->insert("Coral" , "#FF7F50");
  this->insert("DarkOrange" , "#FF8C00");
  this->insert("LightSalmon" , "#FFA07A");
  this->insert("Orange" , "#FFA500");
  this->insert("LightPink" , "#FFB6C1");
  this->insert("Pink" , "#FFC0CB");
  this->insert("Gold" , "#FFD700");
  this->insert("PeachPuff" , "#FFDAB9");
  this->insert("NavajoWhite" , "#FFDEAD");
  this->insert("Moccasin" , "#FFE4B5");
  this->insert("Bisque" , "#FFE4C4");
  this->insert("MistyRose" , "#FFE4E1");
  this->insert("BlanchedAlmond" , "#FFEBCD");
  this->insert("PapayaWhip" , "#FFEFD5");
  this->insert("LavenderBlush" , "#FFF0F5");
  this->insert("SeaShell" , "#FFF5EE");
  this->insert("Cornsilk" , "#FFF8DC");
  this->insert("LemonChiffon" , "#FFFACD");
  this->insert("FloralWhite" , "#FFFAF0");
  this->insert("Snow" , "#FFFAFA");
  this->insert("Yellow" , "#FFFF00");
  this->insert("LightYellow" , "#FFFFE0");
  this->insert("Ivory" , "#FFFFF0");
  this->insert("White" , "#FFFFFF");
  //----- Insert some Visual Distinct Colors
  InsertVisualDistinctColor("Red");
  InsertVisualDistinctColor("Blue");
  InsertVisualDistinctColor("Magenta");
  InsertVisualDistinctColor("Orange");
  InsertVisualDistinctColor("Green");
  InsertVisualDistinctColor("Pink");
  InsertVisualDistinctColor("Brown");
  InsertVisualDistinctColor("Purple");
  InsertVisualDistinctColor("Yellow");
  InsertVisualDistinctColor("Cyan");
  InsertVisualDistinctColor("Gray");
  InsertVisualDistinctColor("Olive");
  InsertVisualDistinctColor("Beige");
  InsertVisualDistinctColor("Lime");
  InsertVisualDistinctColor("Mint");
  InsertVisualDistinctColor("Teal");
  InsertVisualDistinctColor("Navy");
  InsertVisualDistinctColor("Lavender");
}

ColorClass::ColorClass(void){ this->init();  }


// ColorClass::ColorClass(string filename)
// {
//   ifstream file;
//   string line,colorname,rgb;
//   colorelement e;
//   this->MaxColors = 500;
//   this->vecolor.reserve(this->MaxColors);
//   this->init();
//   file.open(filename);
//   if (!file){cout<<"Error: Could not open file \""<<filename<<"\".\n";exit(0);}
//   while(std::getline(file, line)) {
//     if (this->vecolor.size() >= this->MaxColors) // In the future, change to resize
//       {cout << "Error: Cannot insert more colors. Redefine MaxColors.\n"; exit(0);}
//     if (is_comment(line,"#")) continue;
//     replace(line.begin(),line.end(), '\t',' ');
//     istringstream token(line); // Declare an input string stream.
//     getline(token, colorname,' ');
//     getline(token, rgb,' ');
//     this->insert(colorname,rgb);
//   }
//   this->vecolor.resize(vecolor.size());
//   sort (this->vecolor.begin(), this->vecolor.end(), compare_color);
//   for (int i=0;i<this->vecolor.size();i++)this->code[this->vecolor[i].colorname]=i;
// }

void ColorClass::list()
{
  std::map<string,int>::iterator it = this->code.begin();
  cout << "List of colors" << endl;
  while (it!=code.end()){cout<<
      setw(30)<<it->first<<
      setw(5)<<it->second<<
      setw(10)<<ColorRGB(it->second)<<
      endl;it++;}
}

void ColorClass::print()
{
  std::cout << "List of colors" << endl <<
    setw(5) << "Code" << " " <<
    setw(30) << "Color_Name" <<
    setw(10) << "RGB" << endl;
  for (int i=0;i<this->vecolor.size();i++) {
    std::cout << setw(5) << i << " " <<
      setw(30) << this->vecolor[i].colorname <<
      setw(10) << this->vecolor[i].rgb <<  endl;}
}
