// =============================================================
//
//  Some useful general routines
//
//  Send comments/corrections to Flavio K. Miyazawa.
//
// =============================================================
#include <iostream>
#include "myutils.h"
#include <sys/stat.h>
#include <iomanip>
using namespace std;
#include <sstream>
#include <cstring>

// The code below is divided into sections:
//     * Global Variables
//     * Dealing with PDF files
//     * Functions to test values
//     * Type Conversion Routines
//     * Utility functions
//     * Routines to read a table from text file
//     *
//     *

// Some global variables
double MY_EPS;
double MY_INF;
myutilsdefaultparameters myutilsparameters;

//====================================================================================
//     * Dealing with PDF files
// Set the program to view pdf files
void set_pdfreader(string programname) { myutilsparameters.pdfreader = programname; }

// To see a pdf file. It uses the pdf reader defined by set_pdfreader.
int view_pdffile(string filename)
{
  char cmd[MAXCOMMANDSIZE];
  sprintf(cmd, "%s %s", myutilsparameters.pdfreader.c_str(), filename.c_str());
  system(cmd);
  return (0);
}

bool VectorIsInteger(vector<double> &v)
{
  for (vector<double>::iterator it = v.begin(); it != v.end(); ++it)
    if (IsFrac(*it))
      return false;
  return (true);
}

//===================================================================================
//    ROUTINES FOR TIME MANIPULATION

/* return time in seconds since 1/jan/1970*/
long time70(void)
{
  return ((long)time(0));
}

// given a time in seconds, print corresponding time in days, hours, minutes and seconds
void printtime(long t)
{
  long dias, horas, minutos, segundos;
  if (t < 1)
  {
    printf("<1 second");
    return;
  }
  minutos = (long)t / 60;
  segundos = t - minutos * 60;
  horas = (long)minutos / 60;
  minutos = minutos - horas * 60;
  dias = (long)horas / 24;
  horas = horas - dias * 24;
  if (dias == 1)
    printf("%ld day", dias);
  if (dias > 1)
    printf("%ld days", dias);
  if (dias && (horas || minutos || segundos))
    printf(", ");
  if (horas == 1)
    printf("%ld hour", horas);
  if (horas > 1)
    printf("%ld hours", horas);
  if (horas && (minutos || segundos))
    printf(", ");
  if (minutos == 1)
    printf("%ld minute", minutos);
  if (minutos > 1)
    printf("%ld minutes", minutos);
  if (minutos && segundos)
    printf(", ");
  if (segundos == 1)
    printf("%ld second", segundos);
  if (segundos > 1)
    printf("%ld seconds", segundos);
}

// given a time in seconds, generate a string containing the time in days, hours, minutes and seconds
void sprinttime(char *s, long t)
{
  long dias, horas, minutos, segundos;
  char aux[100];
  if (t < 1)
  {
    sprintf(s, "%s", "<1 second");
    return;
  }
  minutos = (long)t / 60;
  segundos = t - minutos * 60;
  horas = (long)minutos / 60;
  minutos = minutos - horas * 60;
  dias = (long)horas / 24;
  horas = horas - dias * 24;
  s[0] = '\0';
  if (dias == 1)
  {
    sprintf(aux, "%ld day", dias);
    strcat(s, aux);
  }
  if (dias > 1)
  {
    sprintf(aux, "%ld days", dias);
    strcat(s, aux);
  }
  if (dias && (horas || minutos || segundos))
    strcat(s, ", ");
  if (horas == 1)
  {
    sprintf(aux, "%ld hour", horas);
    strcat(s, aux);
  }
  if (horas > 1)
  {
    sprintf(aux, "%ld hours", horas);
    strcat(s, aux);
  }
  if (horas && (minutos || segundos))
    strcat(s, ", ");
  if (minutos == 1)
  {
    sprintf(aux, "%ld minute", minutos);
    strcat(s, aux);
  }
  if (minutos > 1)
  {
    sprintf(aux, "%ld minutes", minutos);
    strcat(s, aux);
  }
  if (minutos && segundos)
    strcat(s, ", ");
  if (segundos == 1)
  {
    sprintf(aux, "%ld second", segundos);
    strcat(s, aux);
  }
  if (segundos > 1)
  {
    sprintf(aux, "%ld seconds", segundos);
    strcat(s, aux);
  }
}

// print the time in the format dd,hh:mm:ss
void shortprinttime(long t)
{
  long dias, horas, minutos, segundos;
  minutos = (long)t / 60;
  segundos = t - minutos * 60;
  horas = (long)minutos / 60;
  minutos = minutos - horas * 60;
  dias = (long)horas / 24;
  horas = horas - dias * 24;
  if (dias > 0)
    printf("%1ldd,%2ldh:%2ldm:%2lds", dias, horas, minutos, segundos);
  else if (horas > 0)
    printf("%1ldh:%2ldm:%2lds", horas, minutos, segundos);
  else if (minutos > 0)
    printf("%2ldm:%2lds", minutos, segundos);
  else
    printf("%2lds", segundos);
}
// ========================================================================
//     * Routines to read a table of strings from a text file
//
// read a table containing a header (first line that is not a comment) and then,
// read nr rows, according to the header.
//
StringTable::StringTable(int nr, ifstream &file)
{
  string word, h;
  this->nrows = nr;
  this->line.reserve(nr);
  line.resize(nr);
  this->ncols = 0;
  while ((std::getline(file, h)) && is_comment(h, "#"))
    ;
  replace(h.begin(), h.end(), '\t', ' ');
  istringstream token(h); // Declare an input string stream.
  // read the header
  while (getline(token, word, ' '))
  {
    trim(word);
    if (word.length() == 0)
      continue;
    lowercase(word);
    this->header.push_back(word);
    this->columnsize.push_back(word.length());
    this->ncols++;
  }
  // read the rows
  for (int i = 0; i < nrows; i++)
  {
    while ((std::getline(file, h)) && is_comment(h, "#"))
      ;
    replace(h.begin(), h.end(), '\t', ' ');
    istringstream token2(h); // Declare an input string stream.
    int c = 0;
    this->line[i].reserve(this->ncols);
    this->line[i].resize(this->ncols);
    while (getline(token2, word, ' '))
    {
      trim(word);
      if (word.length() == 0)
        continue;
      if (c == this->ncols)
      {
        cout << i + 1 << ": \"" << h << "\"" << endl;
        cout << "Number of cols, in line " << i + 1 << " is larger than in header." << endl;
        exit(1);
      }
      this->line[i][c] = word;
      if ((int) word.length() > this->columnsize[c])
        this->columnsize[c] = word.length();
      c++;
    }
    if (c < this->ncols)
    {
      cout << "Number of cols, in line " << i + 1 << " is smaller than in header." << endl;
      //exit(1);
    }
  }
  return;
}
StringTable::StringTable(int _nrows, int _ncols)
{
  string word, h;
  this->nrows = _nrows;
  this->ncols = _ncols;
  this->header.reserve(_ncols);
  this->header.resize(_ncols);
  this->line.reserve(_nrows);
  line.resize(_nrows);
  for (int i = 0; i < _nrows; i++)
  {
    this->line[i].reserve(_ncols);
    this->line[i].resize(_ncols);
  }
  return;
}

string StringTable::first(string columnname)
{
  int col = this->columnindex(columnname);
  if (col != -1)
    return (line[0][col]);
  else
  {
    cout << "Table does not contain column name \"" << columnname << "\"." << endl;
    exit(1);
  }
}

int StringTable::firstint(string columnname)
{
  int col = this->columnindex(columnname);
  if (col != -1)
    return (StringToInt(line[0][col]));
  else
  {
    cout << "Table does not contain column name \"" << columnname << "\"." << endl;
    exit(1);
  }
}
double StringTable::firstdouble(string columnname)
{
  int col = this->columnindex(columnname);
  if (col != -1)
    return (StringToDouble(line[0][col]));
  else
  {
    cout << "Table does not contain column name \"" << columnname << "\"." << endl;
    exit(1);
  }
}

bool StringTable::readcolumn(string columnname, vector<string> &scol)
{
  int col = this->columnindex(columnname);
  if (col == -1)
    return (false);
  scol.reserve(this->nrows);
  scol.resize(this->nrows);
  for (int i = 0; i < this->nrows; i++)
    scol[i] = this->line[i][col];
  return (true);
}

bool StringTable::readcolumn(string columnname, vector<int> &scol)
{
  int col = this->columnindex(columnname);
  if (col == -1)
    return (false);
  scol.reserve(this->nrows);
  scol.resize(this->nrows);
  for (int i = 0; i < this->nrows; i++)
    scol[i] = StringToInt(this->line[i][col]);
  return (true);
}

bool StringTable::readcolumn(string columnname, vector<double> &scol)
{
  int col = this->columnindex(columnname);
  if (col == -1)
    return (false);
  scol.reserve(this->nrows);
  scol.resize(this->nrows);
  for (int i = 0; i < this->nrows; i++)
    scol[i] = StringToDouble(this->line[i][col]);
  return (true);
}

void StringTable::print(void)
{
  for (int j = 0; j < this->ncols; j++)
  {
    for (int k = 0; k < this->columnsize[j]; k++)
      cout << "-";
    cout << "-+ ";
  }
  cout << endl;
  for (int j = 0; j < this->ncols; j++)
    cout << std::setw(this->columnsize[j]) << header[j] << " | ";
  cout << endl;
  for (int j = 0; j < this->ncols; j++)
  {
    for (int k = 0; k < this->columnsize[j]; k++)
      cout << "-";
    cout << " + ";
  }
  cout << endl;
  for (int i = 0; i < this->nrows; i++)
  {
    for (int j = 0; j < this->ncols; j++)
      cout << std::setw(this->columnsize[j]) << this->line[i][j] << " | ";
    cout << endl;
  }
  for (int j = 0; j < this->ncols; j++)
  {
    for (int k = 0; k < this->columnsize[j]; k++)
      cout << "-";
    cout << " + ";
  }
  cout << endl;
}
