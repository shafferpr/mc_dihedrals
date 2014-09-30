#include <Two_d_grid.h> //includes the 2d grid class for handling everything you want to do with a 2d grid
#include <random>
#include <chrono>

using namespace std;
int Nbiases=0;
int gridx=200;
int gridy=200;
int N_cvs=0;
double step_size=0.1;
double beta=5;
double energy=0;
int Nsteps=1000;
int Nsweeps=10000;
vector <Two_d_grid> bias_grid;

void initialize(string);
void run_mc();

int main(int argc, char *argv[]){
  string biasfilename;
  if(argc>1){
    N_cvs = atoi( argv[1]);
    Nbiases = atoi (argv[2]);
    biasfilename= argv[3];
  }
  cout << "hello\n";
  initialize(biasfilename);
  run_mc();
  
  return 0; 
}

void run_mc(){
  vector <double> positions(20,0.0);
  double pos_temp=0;
  double energy_delta=0;
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  uniform_real_distribution<double> real_distribution(0.0,1.0);
  uniform_int_distribution<int> distribution(0,N_cvs-1);
  int dihedral_label=0;
  for(int i=0; i<=Nsweeps; i++){
    for(int j=0; j<=Nsteps; j++){
      dihedral_label = distribution(generator);
      pos_temp = positions[dihedral_label] + real_distribution(generator)*step_size;
      if(dihedral_label>0 && dihedral_label< N_cvs-1){
	energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(pos_temp, positions[dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(positions[dihedral_label],positions[dihedral_label+1]);
	energy_delta = bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1], pos_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1],positions[dihedral_label]);
      }
      else if(dihedral_label=0){	
	energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(pos_temp, positions[dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(positions[dihedral_label],positions[dihedral_label+1]);
      }
      else{
	energy_delta = bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1], pos_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1],positions[dihedral_label]);
	if(energy_delta<0 || real_distribution(generator)<=exp(-beta*energy_delta)){
	  energy += energy_delta;
	  positions[dihedral_label]=pos_temp;
	}
	else{
	  pos_temp=pos_temp;
	}
      }
    }
  }

}


void initialize(string biasfname){
  cout << "hello\n";
  Two_d_grid gridtemplate(gridx, gridy);
//Nbiases is the number of 2d_grids you want to make, gridx, and gridy are the number of gridpoints you are using
  cout << "hello\n";
  for(int i=0; i<Nbiases; i++){
    cout << i << "a\n";
    bias_grid.push_back(gridtemplate);
    cout << i << "b\n";
    bias_grid[i].readfiles(i,biasfname);
    cout << i << "c\n";
    bias_grid[i].setindices();
    cout << i << "d\n";
  }

}
