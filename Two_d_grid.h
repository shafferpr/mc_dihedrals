#ifndef __TWO_D_GRID_H_INCLUDED__
#define __TWO_D_GRID_H_INCLUDED__

#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;


class Two_d_grid
{
  double xrange, yrange;
  double xmin, ymin;
  double xspacing, yspacing;
 public:
  std::vector< vector <double> > GridValuesAndCoordinates;
  std::vector< vector <double> > GridValuesByIndex;
  std::vector<int> grid_dimensions;
  int grid_size;
  void readfiles(int,string);
  void setindices();
  double getvalue_linearinterpolation(double, double);
  Two_d_grid(int,int);
  Two_d_grid();
};

/*
Two_d_grid::Two_d_grid()
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
  }*/

Two_d_grid::Two_d_grid(int a, int b)
{
  grid_dimensions.resize(2);
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

void Two_d_grid::readfiles(int pair_identifier, string biasfilename_)
{
  string s;
  string line;
  int counter=0;
  double dum1, dum2=0;
  stringstream ss;
  vector <double> values(3);
  char values_char[10];
  s=std::to_string(pair_identifier);
  ifstream bias_file (biasfilename_+s, ifstream::in);
  //bias_file.open(biasfilename_+s);
  while( getline(bias_file, line)){
    if(!line.empty()){
      ss << line;
      //ss.str(line);
      //cout << ss << "\n";
      /*for(int i=0; i<3; i++){
	ss.getline(values_char);
	//cout << values_char[0] << "\n";
	values[i]=atof(values_char);
	}*/
      ss >> values[0] >> values[1] >> values[2] >> dum1 >> dum2;
      values[2]=-values[2];
      counter += 1;
      ss.clear();
      GridValuesAndCoordinates.push_back(values);
      //if(pair_identifier==3 && counter>39990)
      //cout << values[0]<< " " << values[1] << " " << values[2] << "\n";
    }
  }
  bias_file.close();
}

void Two_d_grid::setindices()
{
  for(int i=0; i<grid_size; i++){
    double x= GridValuesAndCoordinates[i][0]-xmin;
    int index1=(x+0.0001)/xspacing;
    double y=GridValuesAndCoordinates[i][1]-ymin;
    int index2=(y+0.0001)/yspacing;
    GridValuesByIndex[index1][index2]=GridValuesAndCoordinates[i][2];
    //if(index1==100 && index2==101)
    //cout << GridValuesAndCoordinates[i][0] <<" " << GridValuesAndCoordinates[i][1] << " " << GridValuesByIndex[index1][index2] <<" " << xspacing << "\n";
  }

}

double Two_d_grid::getvalue_linearinterpolation(double x, double y)
{
  signed int index1 = (x-xmin+0.0001)/xspacing;
  signed int index2 = (y-ymin+0.0001)/yspacing;
  double x1 = xmin+xspacing*index1;
  double y1 = ymin+yspacing*index2;
  vector <double> point1(3);
  vector <double> point2(3);
  vector <double> point3(3);
  vector <double> cross_product(3);
  vector <double> point2_point1(3);
  vector <double> point3_point1(3);
  //if(index1>198 || index1 <0 || index2>199 || index2<0)
  // printf("%d %d\n", index1, index2);
  point1[0]=x1; point1[1]=y1; point1[2]=GridValuesByIndex[index1][index2];
  point2[0]=x1+xspacing; point2[1]=y1;
  point3[0]=x1; point3[1]=y1+yspacing;

  if(index1==(grid_dimensions[0]-1) && index2==(grid_dimensions[1]-1)){
    point2[2]=GridValuesByIndex[0][index2];
    point3[2]=GridValuesByIndex[index1][0];
  }
  else if(index1==(grid_dimensions[0]-1)){
    point2[2]=GridValuesByIndex[0][index2];
  }
  else if(index2==(grid_dimensions[1]-1)){
    point3[2]=GridValuesByIndex[index1][0];
  }
  else{
    point2[2]=GridValuesByIndex[index1+1][index2];
    point3[2]=GridValuesByIndex[index1][index2+1]; 
  }
  //printf("a %f b %f c %f %d %d\n", point1[0], point1[1], point1[2], index1, index2);
  //printf("a %f b %f c %f\n", point2[0], point2[1], point2[2]);
  //printf("a %f b %f c %f\n", point3[0], point3[1], point3[2]);
  point2_point1[0]=point2[0]-point1[0];point2_point1[1]=point2[1]-point1[1];point2_point1[2]=point2[2]-point1[2];
  point3_point1[0]=point3[0]-point1[0];point3_point1[1]=point3[1]-point1[1];point3_point1[2]=point3[2]-point1[2];
  cross_product[0]=point2_point1[1]*point3_point1[2]-point2_point1[2]*point3_point1[1];
  cross_product[1]=point2_point1[2]*point3_point1[0]-point2_point1[0]*point3_point1[2];
  cross_product[2]=point2_point1[0]*point3_point1[1]-point2_point1[1]*point3_point1[0];
  double r2= sqrt(pow(cross_product[0],2)+pow(cross_product[1],2)+pow(cross_product[2],2));
  cross_product[0] /= r2; cross_product[1] /= r2; cross_product[2] /= r2;
  double d= -(cross_product[0]*point1[0]+cross_product[1]*point1[1]+cross_product[2]*point1[2]);
  double z = (-d - x*cross_product[0] - y*cross_product[1])/cross_product[2];
  //printf("z %f\n", z);
  //printf("%f %f %f\n", d/cross_product[2], cross_product[0]/cross_product[2], cross_product[1]/cross_product[2]);
  return z;
}

#endif
