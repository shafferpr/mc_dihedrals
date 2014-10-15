#include <Two_d_grid.h> //includes the 2d grid class for handling everything you want to do with a 2d grid
#include <random>
#include <chrono>
#include <armadillo>

using namespace arma;
using namespace std;
int Nbiases=0;
int gridx=200;
int gridy=200;
int N_cvs=0;
double step_size=0.1;
double beta=0.401606;
double energy=0;
int Nbackbone=0;
int Nsteps=5500;
int Nsweeps=3000000;
//int Nsteps=1;
//int Nsweeps=1;
vector <Two_d_grid> bias_grid;

void initialize(string);
void run_mc();
void print(int, vector<double> &, double);
void backbone(int, vector<double> &);

int main(int argc, char *argv[]){
  string biasfilename;
  if(argc>1){
    N_cvs = atoi( argv[1]);
    Nbiases = atoi (argv[2]);
    biasfilename= argv[3];
  }
  //cout << "hello\n";
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

      if(pos_temp>xmax-0.01){
	pos_temp=xmin+(pos_temp-xmax);
      }
      else if(pos_temp<xmin){
	pos_temp=xmax-(xmin-pos_temp)-0.001;
      }

      if(dihedral_label>0 && dihedral_label< N_cvs-1){
	energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(pos_temp, positions[dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(positions[dihedral_label],positions[dihedral_label+1]);

	energy_delta += bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1], pos_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1],positions[dihedral_label]);

      }
      else if(dihedral_label==0){	
	energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(pos_temp, positions[dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(positions[dihedral_label],positions[dihedral_label+1]);

      }
      else{
	energy_delta = bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1], pos_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[dihedral_label-1],positions[dihedral_label]);

      }

      if(energy_delta<0 || real_distribution(generator)<=exp(-beta*energy_delta)){
	energy += energy_delta;
	positions[dihedral_label]=pos_temp;

      }
      else{
	//cout << "rejected\n";
	pos_temp=pos_temp;
      }
    }
    backbone(i, positions);  
  }

}

void backbone(int counter, vector<double> & allpositions){
  vector< vector <double> > backbone;
  vector< double > bond_distances(Nbackbone+6);
  vector< double > bond_angles(Nbackbone+6);
  vector< double > dihedral_angles(Nbackbone+6);
  ofstream backbonefile;
  string bbstring ("backbone_atoms.xyz");
  backbonefile.open(bbstring,ios::out);
  backbone.resize(Nbackbone);

  /*for(int k=0; k<Nbackbone; k++){
    backbone[k].resize(3);
  }

  backbone[0][0]=0;backbone[0][1]=0;backbone[0][2]=0;
  bond_distances[0]=0.133386;
  for(int k=0; k<=N_cvs/2; k++){
    bond_distances[1+3*k]=0.147101;
    bond_distances[2+3*k]=0.1546;
    bond_distances[3+3*k]=0.133649;
  }
  bond_angles[0]=2.1957;
  for(int k=0; k<=N_cvs/2; k++){
    bond_angles[1+3*k]=2.0031;
    bond_angles[2+3*k]=2.0653;
    bond_angles[3+3*k]=2.17011;
  }
  int j1=0;
  allpositions[0]=-1.5;
  allpositions[1]=0.01713;
  allpositions[2]=0.8095;
  allpositions[3]=-1.6332;
  allpositions[4]=1.642;
  allpositions[5]=-1.585;
  allpositions[6]=2.04;
  allpositions[7]=-1.4;
  allpositions[8]=3.0;
  allpositions[9]=-0.5;

  allpositions[0]=3.1415926;
  allpositions[1]=3.1415926;
  allpositions[2]=3.1415926;
  allpositions[3]=3.1415926;
  allpositions[4]=3.1415926;
  allpositions[5]=3.1415926;
  allpositions[6]=3.1415926;
  allpositions[7]=3.1415926;
  allpositions[8]=3.1415926;
  allpositions[9]=3.1415926;

  for(int k=0; k<Nbackbone/3; k++){
    dihedral_angles[3*k+0]=allpositions[j1];j1++;
    cout <<dihedral_angles[3*k+0] << "\n";
    dihedral_angles[3*k+1]=allpositions[j1];j1++;
    cout <<dihedral_angles[3*k+1] << "\n";
    dihedral_angles[3*k+2]=3.14159265;
    cout <<dihedral_angles[3*k+2] << "\n";
  }
  vec terminal_vector(4);
  mat transformation_matrix(4,4,fill::eye);
  mat new_matrix(4,4,fill::eye);
  new_matrix(0,0)=-1;
  new_matrix(2,2)=-1;
  new_matrix(0,3)=-bond_distances[0];
  terminal_vector(0)=0; terminal_vector(1)=0; terminal_vector(2)=0; terminal_vector(3)=1;
  
  backbone[1][0]=-bond_distances[0];backbone[1][1]=0; backbone[1][2]=0;
  
  transformation_matrix=transformation_matrix*new_matrix;
  new_matrix(0,0) = -cos(bond_angles[0]); new_matrix(0,1) = -sin(bond_angles[0]); new_matrix(0,2)=0; new_matrix(0,3)= -bond_distances[1]*cos(bond_angles[0]);
  new_matrix(1,0) = sin(bond_angles[0]); new_matrix(1,1) = -cos(bond_angles[0]); new_matrix(1,2)=0; new_matrix(1,3)= bond_distances[1]*sin(bond_angles[0]);
  new_matrix(2,2)=1;

  transformation_matrix=transformation_matrix*new_matrix;
  backbone[2][0]=transformation_matrix(0,3);
  backbone[2][1]=transformation_matrix(1,3);
  backbone[2][2]=transformation_matrix(2,3);

  backbonefile <<"C "<<backbone[0][0] <<" "<<backbone[0][1]<<" "<<backbone[0][2]<<"\n";
  backbonefile <<"N "<<backbone[1][0] <<" "<<backbone[1][1]<<" "<<backbone[1][2]<<"\n";
  backbonefile <<"C "<<backbone[2][0] <<" "<<backbone[2][1]<<" "<<backbone[2][2]<<"\n";

  double cos_ba=0; double sin_ba=0;
  double cos_da=0; double sin_da=0;
  for(int k=3; k<Nbackbone; k++){
    cos_ba=cos(bond_angles[k-2]); sin_ba=sin(bond_angles[k-2]);
    cos_da=cos(dihedral_angles[k-3]); sin_da=sin(dihedral_angles[k-3]);
    
    new_matrix(0,0)=-cos_ba; new_matrix(0,1)=-sin_ba; new_matrix(0,2)=0, new_matrix(0,3)= -bond_distances[k-1]*cos_ba;
    new_matrix(1,0)=sin_ba*cos_da; new_matrix(1,1)=-cos_ba*cos_da; new_matrix(1,2)= -sin_da; new_matrix(1,3) = bond_distances[k-1]*sin_ba*cos_da;
    new_matrix(2,0)=sin_ba*sin_da; new_matrix(2,1)=-cos_ba*sin_da; new_matrix(2,2)= cos_da; new_matrix(2,3) = bond_distances[k-1]*sin_ba*sin_da;
    new_matrix(3,0)=0; new_matrix(3,1)=0; new_matrix(3,2)=0; new_matrix(3,3)=1;
    transformation_matrix=transformation_matrix*new_matrix;

    backbone[k][0]=transformation_matrix(0,3);
    backbone[k][1]=transformation_matrix(1,3);
    backbone[k][2]=transformation_matrix(2,3);

    if(k%3==1){
      backbonefile <<"N " <<backbone[k][0] <<" "<<backbone[k][1]<<" "<<backbone[k][2]<<"\n";
    }
    else{
      backbonefile <<"C " <<backbone[k][0] <<" "<<backbone[k][1]<<" "<<backbone[k][2]<<"\n";
    }
    
  }

  double rx, ry, rz=0;
  //rx=backbone[14][0]-backbone[2][0];
  //ry=backbone[14][1]-backbone[2][1];
  //rz=backbone[14][2]-backbone[2][2];
  double r2 = sqrt(rx*rx+ry*ry+rz*rz);*/
  double r2=0;

  print(counter, allpositions, r2);
  //cout << "hello 7\n";
  //cout << r2 <<" "<< bond_distances[3]<<" r2\n";
}


void initialize(string biasfname){
  Two_d_grid gridtemplate(gridx, gridy);
//Nbiases is the number of 2d_grids you want to make, gridx, and gridy are the number of gridpoints you are using
  for(int i=0; i<Nbiases; i++){
    bias_grid.push_back(gridtemplate);
    bias_grid[i].readfiles(i,biasfname);
    bias_grid[i].setindices();
  }
  cout <<  bias_grid[0].GridValuesByIndex[20][20] <<"hello \n";
  cout << bias_grid[0].GridValuesAndCoordinates[0][2] << "\n";
  Nbackbone=3*N_cvs/2 +1;
}

void print(int counter, vector<double> & allpositions, double r2){
  ofstream colvarfile;
  colvarfile.open("colvar_nopt.data", ios_base::app);

  for(int i=0; i<N_cvs; i++){
    colvarfile << allpositions[i];
    colvarfile << " ";
  }

  colvarfile << r2;
  colvarfile << "\n";

  colvarfile.close();

}
