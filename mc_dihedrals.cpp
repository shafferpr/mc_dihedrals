#include <Two_d_grid.h> //includes the 2d grid class for handling everything you want to do with a 2d grid
#include <random>
#include <chrono>

using namespace std;
int Nbiases=0;
int gridx=200;
int gridy=200;
int N_cvs=0;
double step_size=0.1;
double beta=0.401606;
double energy=0;
int Nsteps=1000;
int Nsweeps=30000;
vector <Two_d_grid> bias_grid;

void initialize(string);
void run_mc();
void print(int, vector<double> &);

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
  double xmax=3.141592654;
  double xmin=-3.14159264;
  for(int i=0; i<=Nsweeps; i++){
    for(int j=0; j<=Nsteps; j++){
      dihedral_label = distribution(generator);
      pos_temp = positions[dihedral_label] + (real_distribution(generator)-0.5)*step_size;
      //pos_temp=3.14159;
      //cout << pos_temp << "\n";
      if(pos_temp>xmax-0.01){
	pos_temp=xmin+(pos_temp-xmax);
      }
      else if(pos_temp<xmin){
	pos_temp=xmax-(xmin-pos_temp)-0.001;
      }
      //cout << j << "\n";
      if(dihedral_label>0 && dihedral_label< N_cvs-1){
	energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(pos_temp, positions[dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(positions[dihedral_label],positions[dihedral_label+1]);
	energy_delta = bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1], pos_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1],positions[dihedral_label]);
	//cout << energy_delta << " in between  " << "\n";
      }
      else if(dihedral_label==0){	
	energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(pos_temp, positions[dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(positions[dihedral_label],positions[dihedral_label+1]);
	//cout << energy_delta <<" "<<bias_grid[dihedral_label].getvalue_linearinterpolation(positions[dihedral_label],positions[dihedral_label+1]) << " 0  " << "\n";
      }
      else{
	energy_delta = bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1], pos_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1],positions[dihedral_label]);
	//cout << energy_delta <<" ncv " << "\n";
      }
      if(energy_delta<0 || real_distribution(generator)<=exp(-beta*energy_delta)){
	energy += energy_delta;
	positions[dihedral_label]=pos_temp;
	//cout << energy_delta <<" accepted " << "\n";
      }
      else{
	//cout << "rejected\n";
	pos_temp=pos_temp;
      }
    }
    print(i, positions);
  }

}


void initialize(string biasfname){
  Two_d_grid gridtemplate(gridx, gridy);
//Nbiases is the number of 2d_grids you want to make, gridx, and gridy are the number of gridpoints you are using
  for(int i=0; i<Nbiases; i++){
    bias_grid.push_back(gridtemplate);
    bias_grid[i].readfiles(i,biasfname);
    bias_grid[i].setindices();
  }

}

void print(int counter, vector<double> & allpositions){
  ofstream colvarfile;
  colvarfile.open("colvar.data", ios_base::app);
  for(int i=0; i<N_cvs; i++){
    colvarfile << allpositions[i];
    colvarfile << " ";
  }
  colvarfile << "\n";
  colvarfile.close();
}
