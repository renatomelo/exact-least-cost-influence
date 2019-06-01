// =============================================================
//
//  Some useful general routines
//
//  Send comments/corrections to Flavio K. Miyazawa.
//
// =============================================================
#ifndef MYUTILS_DEFINE
#define MYUTILS_DEFINE
#include <string>
#include <vector>
#include <time.h>
#include <math.h>
using namespace std;
//#include <ostream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

//================================================================================================================
// SOME CONSTANTS
//
// maximum number of characters used by a command used by the system routine.
#define MAXCOMMANDSIZE 1000
// used to verify if number are the same, within some error of MY_EPS

//================================================================================================================
//    ROUTINES TO DEAL WITH PDF

//==========================================================================================================
//     * Global Variables
//
class myutilsdefaultparameters
{
public:
  myutilsdefaultparameters();
  string pdfreader;
  double *P_MY_EPS;
  double *P_MY_INF;
};
extern double MY_EPS;
extern double MY_INF;

// Constructor for the global variables declared in myutils.
// Put all here, so we can have visibility of all global variables
inline myutilsdefaultparameters::myutilsdefaultparameters()
{
  pdfreader = "evince"; // Program to open a pdf file in Mac OS X
  // You can change the term "evince" by any other pdf viewer.
  // But the program must run in a terminal with the sintax
  // <pdfreader>  <pdffile>
  // Example:
  //             evince  file.pdf
  //             open    file.pdf
  //             xpdf    file.pdf
  P_MY_EPS = &MY_EPS;
  P_MY_INF = &MY_INF;
  *P_MY_EPS = 0.00001;
  *P_MY_INF = 1000000000000.0;
}
extern myutilsdefaultparameters myutilsparameters;

inline void set_MY_EPS(double neweps) { *(myutilsparameters.P_MY_EPS) = neweps; }
inline void set_MY_INF(double newinf) { *(myutilsparameters.P_MY_INF) = newinf; }

// set default pdf reader used in some programs
// For example, if it is possible to open a pdf file
// with the command "open file.pdf"
// then, you can use the following command to define
// the "open" as a pdf reader
// set_pdfreader("open");
void set_pdfreader(string programname);

// To see a pdf file. It uses the pdf reader defined by set_pdfreader.
int view_pdffile(string filename);

// Only to see if a file exists. It (tries to) open and close the file.
//
// bool FileExists(const char *filename)
// {
//   bool r;
//   ifstream ifile(filename);
//   r = ifile;
//   ifile.close();
//   return r;
// }
//

inline bool is_space(char c)
{
  if ((c == ' ') || (c == '\t'))
    return (true);
  return (false);
}
inline bool is_comment(string line)
{
  unsigned int i = 0;
  cout << "line in function: " << line << endl;
  while (i < (line.length()) && (is_space(line[i])))
    i++;
  if (i == line.length())
    return (true); // consider empty line as a comment
  if (line[i] == '#')
    return (true); // first non-null character is #
  return (false);
}
// comment if prefix appear after blanks -- ROUTINE NOT TESTED
inline bool is_comment(string line, string prefix)
{
  unsigned int i = 0;
  while (i < (line.length()) && (is_space(line[i])))
    i++;
  if (i == line.length())
    return (true); // consider empty line as a comment
  if (i + prefix.length() > line.length())
    return (false);
  return (!line.compare(i, prefix.length(), prefix));
}

inline bool IsSuffix(string str, string suf)
{
  return (str.size() >= suf.size() && !str.compare(str.size() - suf.size(), suf.size(), suf));
}

inline bool IsPrefix(string str, string pre)
{
  return (str.size() >= pre.size() && !str.compare(0, pre.size(), pre));
}

inline int hex2int(string s)
{
  int n = s.length(), d = 0;
  for (int i = 0; i < n; i++)
  {
    char c = tolower(s[i]);
    if ((c >= '0') && (c <= '9'))
      d = d * 16 + (c - '0');
    else if ((c >= 'a') && (c <= 'f'))
      d = d * 16 + (c - 'a' + 10);
    else
    {
      cout << "Error: character \"" << s[i] << "\" in string \""
           << s << "\" is not a hexadecimal digit.\n";
      exit(0);
    }
  }
  return (d);
}

inline bool FileExists(const string &filename)
{
  ifstream f(filename.c_str());
  return f.is_open();
}

//====================================================================================================
//     * Type Conversion Routines
// In C++11 we have new convertion functions, but some libraries
// are still not compatible, like the Gurobi Library.

/* float              stof(const string& str, size_t *idx = 0); */
/* double             stod(const string& str, size_t *idx = 0); */
/* long double        stold(const string& str, size_t *idx = 0); */
/* int                stoi(const string& str, size_t *idx = 0, int base = 10); */
/* long               stol(const string& str, size_t *idx = 0, int base = 10); */
/* unsigned long      stoul(const string& str, size_t *idx = 0, int base = 10); */
/* long long          stoll(const string& str, size_t *idx = 0, int base = 10); */
/* unsigned long long stoull(const string& str, size_t *idx = 0, int base = 10); */
/* string to_string(int val); */
/* string to_string(unsigned val); */
/* string to_string(long val); */
/* string to_string(unsigned long val); */
/* string to_string(long long val); */
/* string to_string(unsigned long long val); */
/* string to_string(float val); */
/* string to_string(double val); */
/* string to_string(long double val); */

/* template <typename T> string ToString ( T val ){ */
/*     std::stringstream retval; */
/*     retval << val; */
/*     return retval; */
/* } */

inline char chartolower(char c)
{
  if (c <= 'Z' && c >= 'A')
    return (c - ('Z' - 'z'));
  return (c);
}

inline void lowercase(string &s)
{
  std::transform(s.begin(), s.end(), s.begin(), chartolower);
}

// left trim
static inline string &ltrim(string &str)
{
  str.erase(str.begin(), std::find_if(str.begin(), str.end(),
                                      std::not1(std::ptr_fun<int, int>(std::isspace))));
  return str;
}

// right trim
static inline string &rtrim(string &str)
{
  str.erase(std::find_if(str.rbegin(), str.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace)))
                .base(),
            str.end());
  return str;
}

// left and right trim
static inline string &trim(string &str) { return ltrim(rtrim(str)); }

// convert a double value to string
inline std::string DoubleToString(double val)
{
  std::stringstream out;
  out << val;
  return out.str();
}
//{return to_string(val);} // this is only valid in c++11

// convert a int value to string
inline std::string IntToString(int val)
{
  std::stringstream out;
  out << val;
  return out.str();
}
//{return to_string(val);}

// convert a string to double
inline double StringToDouble(string s)
{
  double d;
  stringstream(s) >> d;
  return (d);
}
//{return stod(s);}

// convert a string to int
inline int StringToInt(string s)
{
  int n;
  stringstream(s) >> n;
  return (n);
}
//{return stoi(s);}

inline void Pause(void)
{
  cout << "Pause";
  std::cin.get();
  cout << "\n";
}

//====================================================================================
//     * Functions to test values
// Return true if variable x is fractional (within a certain small error).
inline bool IsFrac(double x)
{
  double f;
  f = ceil(x - MY_EPS) - x;
  if (f < MY_EPS)
    return (false);
  if (f > 1.0 - MY_EPS)
    return (false);
  return (true); //    eps  <= ceil(x)-x <= 1-eps
}

// Return true if variable x is fractional (within a certain small error).
inline bool IsEqual(double x, double y)
{
  if (x > y + MY_EPS)
    return (false);
  if (x < y - MY_EPS)
    return (false);
  return (true); //    y-eps  <= x <= y+eps
}

// For the next functions, it is supposed that 0 <= X <= 1
inline bool BinaryIsOne(double x) { return ((x > (double)(1.0 - MY_EPS))); }
inline bool BinaryIsZero(double x) { return ((x < MY_EPS)); }
inline bool NonBinary(double x) { return ((x > MY_EPS) && (x < 1.0 - MY_EPS)); }

// verify if a vector contains only integer values
bool VectorIsInteger(vector<double> &v);

//================================================================================================================
//    ROUTINES FOR TIME MANIPULATION

long time70(void);                /* returns the time in seconds since Jan 1, 1970 */
void printtime(long t);           /* print the time t in days, hours, minutes and seconds.*/
void sprinttime(char *s, long t); /* prints the time in the string s Example: 1 hour, 2 minutes, 3 seconds*/
void shortprinttime(long t);      /* prints the time in the string s. Ex.: 11d,22h:33m:44s   */

//==================================================================================
//     * Routines/classes to read a table of strings from a text file
class StringTable
{
public:
  StringTable(int nrows, ifstream &file);
  StringTable(int nrows, int ncols);
  int nrows, ncols;
  vector<string> header;
  vector<int> columnsize;
  vector<vector<string>> line;
  bool readtable(ifstream &file);
  int columnindex(string columnname);    //return column index, or -1 if it does not exists
  string first(string columnname);       // return the first element of the column
  int firstint(string columnname);       // return the first element of the column as int
  double firstdouble(string columnname); // return the first element of the column as double
  bool readcolumn(string columnname, vector<string> &scol);
  bool readcolumn(string columnname, vector<int> &scol);
  bool readcolumn(string columnname, vector<double> &scol);
  bool entry(int row, int col, string &entry);
  bool entry(int row, int col, double &entry);
  bool entry(int row, int col, int &entry);
  void print();
};

//return the column index, or (-1) if it does not exists
inline int StringTable::columnindex(string columnname)
{
  string cname = columnname;
  lowercase(cname);
  for (int col = 0; col < this->ncols; col++)
    if (cname == this->header[col])
      return (col);
  return (-1);
}

inline bool StringTable::entry(int row, int col, string &entry)
{
  if (row < 0)
    return (false);
  if (row >= this->nrows)
    return (false);
  if (col < 0)
    return (false);
  if (col >= this->ncols)
    return (false);
  entry = this->line[row][col];
  return (true);
}

inline bool StringTable::entry(int row, int col, double &entry)
{
  if (row < 0)
    return (false);
  if (row >= this->nrows)
    return (false);
  if (col < 0)
    return (false);
  if (col >= this->ncols)
    return (false);
  entry = StringToDouble(this->line[row][col]);
  return (true);
}

inline bool StringTable::entry(int row, int col, int &entry)
{
  if (row < 0)
    return (false);
  if (row >= this->nrows)
    return (false);
  if (col < 0)
    return (false);
  if (col >= this->ncols)
    return (false);
  entry = StringToInt(this->line[row][col]);
  return (true);
}

// Given a filename that may include its path, this routine
// return the path, ending with '/'
inline string getpath(string filename)
{
  int n = filename.size();
  if (n == 0)
  {
    cout << "Error: Invalid filename." << endl;
    exit(0);
  }
  n--;
  while ((n >= 0) && (filename[n] != '/'))
    n--;
  if (n < 0)
    return ("./");
  else
    return (filename.substr(0, n + 1));
}

#endif
