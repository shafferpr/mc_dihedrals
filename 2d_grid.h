#ifndef __2D_GRID_H_INCLUDED__
#define __2D_GRID_H_INCLUDED__

#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>

using namespace std;


class 2d_grid
{
  double xrange, yrange;
 public:
  vector< vector <double> > GridValuesAndCoordinates;
  vector< vector <double> > GridValuesByIndex;
  vector< int > grid_dimensions;
  int grid_size;
  void readfiles(int,string);
  void setindices();
  2d_grid(int,int);
  2d_grid();
}

2d_grid::2d_grid()
{
  grid_dimensions[0]=100;
  grid_dimensions[1]=100;
  grid_size=grid_dimensions[0]*grid_dimensions[1];
  xrange=2*3.141592654;
  yrange=2*3.141592654;
  GridValuesByIndex.resize(grid_dimensions[0]);
  for(int i=0; i<grid_dimensions[0]; i++){
    GridValuesByIndex[i].resize(grid_dimensions[1]);
  }
}

2d_grid::2d_grid(int a, int b)
{
  grid_dimensions[0]=a;
  grid_dimensions[1]=b;
  grid_size=grid_dimensions[0]*grid_dimensions[1];
  xrange=2*3.141592654;
  yrange=2*3.141592654;
}

void 2d_grid::readfiles(int pair_identifier, string biasfilename)
{
  string s;
  string line;
  stringstream ss;
  vector <double> values(3);
  char values_char[10];
  s=std::to_string(pair_identifier);
  ifstream bias_file;
  bias_file.open(biasfilename_+s);
  while( getline(bias_file, line)){
    if(!line.empty()){
      ss << line;
      for(int i=0; i<2; i++){
	ss.getline(values_char, 4, '  ');
	values[i]=atod(values_char);
      }
      GridValuesAndCoordinates.push_back(values);
    }
  }
}

void 2d_grid::setindices()
{
  double xspacing=xrange/grid_dimensions[0];
  double yspacing=yrange/grid_dimensions[1];
  for(int i=0; i<grid_size; i++){
    double x=GridValuesAndCoordinates[i][0]-xmin;
    int index1=x/spacing;
    double y=GridValuesAndCoordinates[i][1]-ymin;
    int index2=y/spacing;
    GridValuesByIndex[index1][index2]=GridValuesAndCoordinates[i][3];
  }

}
