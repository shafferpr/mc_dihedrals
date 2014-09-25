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
  double xmin, ymin;
  double xspacing, yspacing;
 public:
  vector< vector <double> > GridValuesAndCoordinates;
  vector< vector <double> > GridValuesByIndex;
  vector< int > grid_dimensions;
  int grid_size;
  void readfiles(int,string);
  void setindices();
  double getvalue_linearinterpolation(double, double);
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
  xmin=-3.141592654;
  ymin=-3.141592654;
  xspacing=xrange/grid_dimensions[0];
  yspacing=yrange/grid_dimensions[1];
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
  xmin=-3.141592654;
  ymin=-3.141592654;
  xspacing=xrange/grid_dimensions[0];
  yspacing=yrange/grid_dimensions[1];
  GridValuesByIndex.resize(grid_dimensions[0]);
  for(int i=0; i<grid_dimensions[0]; i++){
    GridValuesByIndex[i].resize(grid_dimensions[1]);
  }
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
  for(int i=0; i<grid_size; i++){
    double x=GridValuesAndCoordinates[i][0]-xmin;
    int index1=x/xspacing;
    double y=GridValuesAndCoordinates[i][1]-ymin;
    int index2=y/yspacing;
    GridValuesByIndex[index1][index2]=GridValuesAndCoordinates[i][3];
    
  }

}

double 2d_grid::getvalue_linearinterpolation(double x, double y)
{
  int index1 = (x-xmin)/xspacing;
  int index2 = (y-ymin)/yspacing;
  double x1 = xmin+xspacing*index1;
  double y1 = ymin+yspacing*index2;
  vector <double> point1(3);
  vector <double> point2(3);
  vector <double> point3(3);
  vector <double> cross_product(3);
  vector <double> point2_point1(3);
  vector <double> point3_point1(3);
  point1[0]=x1; point1[1]=y1; point1[2]=GridValuesByIndex[index1][index2];
  point2[0]=x1+xspacing; point2[1]=y1;
  point3[0]=x1; point3[1]=y1+yspacing;
  if(index1==grid_dimensions[0] && index2==grid_dimensions[1]){
    point2[2]=GridValuesByIndex[0][index2];
    point3[2]=GridValuesByIndex[index1][0];
  }
  else if(index1==grid_dimensions[0]){
    point2[2]=GridValuesByIndex[0][index2];
  }
  else if(index2==grid_dimensions[1]){
    point3[2]=GridValuesByIndex[index1][0];
  }
  else{
    point2[2]=GridValuesByIndex[index1+1][index2];
    point3[2]=GridValuesByIndex[index1][index2+1]; 
  }
  point2_point1[0]=point2[0]-point1[0];point2_point1[1]=point2[1]-point1[1];point2_point1[2]=point2[2]-point1[2];
  point3_point1[0]=point3[0]-point1[0];point3_point1[1]=point3[1]-point1[1];point3_point1[2]=point3[2]-point1[2];
  cross_product[0]=point2_point1[1]*point3_point1[2]-point2_point1[2]*point3_point1[1];
  cross_product[1]=point2_point1[2]*point3_point1[0]-point2_point1[0]*point3_point1[2];
  cross_product[2]=point2_point1[0]*point3_point1[1]-point2_point1[1]*point3_point1[0];
  double r2= sqrt(pow(cross_product[0],2)+pow(cross_product[1],2)+pow(cross_product[2],2));
  cross_product[0] /= r2; cross_product[1] /= r2; cross_product[2] /= r2;
  double d= -(cross_product[0]*point1[0]+cross_product[1]*point1[1]+cross_product[2]*point1[2]);
  double z = (d - x*cross_product[0] - y*cross_product[1])/cross_product[2];
  return z;
}
