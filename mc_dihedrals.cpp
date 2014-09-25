#include <2d_grid.h>

int Nbiases=0;
int gridx=200;
int gridy=200;

int main(){

  void initialize();
  

}

initialize(){
  vector <2d_grid> bias_grid(Nbiases, gridx, gridy);
  for(i=0; i<Nbiases; i++){
    bias_grid[i].readfiles(i,bias);
    bias_grid[i].setindices();
  }

}
