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
vector <double> beta(8);
int Nbackbone=0;
int Nsteps=500;
int Nsweeps=1000000;
vector <Two_d_grid> bias_grid;

void initialize(string);
void run_mc();
void print(int, vector<double> &, double, vector<double> &);
void backbone(int, vector<double> &, vector<double> &);

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
  vector <vector <double> > positions;
  double pos_temp=0;
  double energy_delta=0;
  int Nreplicas=7;
  vector <double> energy(Nreplicas,0);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  uniform_real_distribution<double> real_distribution(0.0,1.0);
  uniform_int_distribution<int> distributionA(0,N_cvs-1);
  uniform_int_distribution<int> distributionB(0,Nreplicas-2);
  int dihedral_label=0;
  int swap_pair=0;
  double xmax=3.141592654;
  double xmin=-3.14159264;
  positions.resize(Nreplicas);
  
  for(int k=0; k<positions.size(); k++){
    positions[k].resize(N_cvs);
  }
  
  for(int i=0; i<=Nsweeps; i++){
    for(int k=0; k<7; k++){
      for(int j=0; j<=Nsteps; j++){
	dihedral_label = distributionA(generator);
	pos_temp = positions[k][dihedral_label] + (real_distribution(generator)-0.5)*step_size;
	if(pos_temp>xmax-0.01){
	  pos_temp=xmin+(pos_temp-xmax);
	}
	else if(pos_temp<xmin){
	  pos_temp=xmax-(xmin-pos_temp)-0.001;
	}
	energy_delta=0;
	if(dihedral_label>0 && dihedral_label< N_cvs-1){
	  energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(pos_temp, positions[k][dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(positions[k][dihedral_label],positions[k][dihedral_label+1]);
	  energy_delta += bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[k][dihedral_label-1], pos_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[k][dihedral_label-1],positions[k][dihedral_label]);

	}
	else if(dihedral_label==0){	
	  energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(pos_temp, positions[k][dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(positions[k][dihedral_label],positions[k][dihedral_label+1]);
	}
	else{
	  energy_delta = bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[k][dihedral_label-1], pos_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[k][dihedral_label-1],positions[k][dihedral_label]);
	}
	if(energy_delta<0 || real_distribution(generator)<=exp(-beta[k]*energy_delta)){
	  energy[k] += energy_delta;
	  positions[k][dihedral_label]=pos_temp;
	  //cout << energy[k]<< "\n";
	}
	else{
	  //cout << "rejected\n";
	  pos_temp=pos_temp;
	}
      }
    }
    swap_pair=distributionB(generator);
    //cout <<swap_pair<< "\n";
    double bf = (beta[swap_pair]-beta[swap_pair+1])*(energy[swap_pair+1]-energy[swap_pair]);
    //cout <<bf << "\n";
    if(bf>0 || real_distribution(generator)<exp(bf)){
      for(int k=0; k<N_cvs; k++){
	double tmp1=positions[swap_pair][k];
	double tmp2=positions[swap_pair+1][k];
	positions[swap_pair+1][k]=tmp1;
	positions[swap_pair][k]=tmp2;
      }
      double tmp1=energy[swap_pair];
      double tmp2=energy[swap_pair+1];
      energy[swap_pair+1]=tmp1;
      energy[swap_pair]=tmp2;
      //if(swap_pair==4)
      //cout << "accepted "<<swap_pair<< "\n";
    }
    if(i%10==0){
      backbone(i, positions[0], energy);
    }
  }
  
}

void backbone(int counter, vector<double> & allpositions, vector<double> & allenergies){
  vector< vector <double> > backbone;
  vector< double > bond_distances(Nbackbone+6);
  vector< double > bond_angles(Nbackbone+6);
  vector< double > dihedral_angles(Nbackbone+6);
  //ofstream backbonefile;
  //string bbstring ("backbone_atoms.xyz");
  //backbonefile.open(bbstring,ios::out);
  backbone.resize(Nbackbone);
  for(int k=0; k<Nbackbone; k++){
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
  /*allpositions[0]=-0.1974;
    allpositions[1]=0.01713;
    allpositions[2]=0.8095;
    allpositions[3]=-1.6332;
    allpositions[4]=1.642;
    allpositions[5]=-1.585;*/
  for(int k=0; k<Nbackbone/3; k++){
    dihedral_angles[3*k+0]=allpositions[j1];j1++;
    dihedral_angles[3*k+1]=allpositions[j1];j1++;
    dihedral_angles[3*k+2]=3.14159265;
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

  //backbonefile <<"C "<<10*backbone[0][0] <<" "<<10*backbone[0][1]<<" "<<10*backbone[0][2]<<"\n";
  //backbonefile <<"N "<<10*backbone[1][0] <<" "<<10*backbone[1][1]<<" "<<10*backbone[1][2]<<"\n";
  //backbonefile <<"C "<<10*backbone[2][0] <<" "<<10*backbone[2][1]<<" "<<10*backbone[2][2]<<"\n";
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

    /*if(k%3==1){
      backbonefile <<"N " <<10*backbone[k][0] <<" "<<10*backbone[k][1]<<" "<<10*backbone[k][2]<<"\n";
    }
    else{
      backbonefile <<"C " <<10*backbone[k][0] <<" "<<10*backbone[k][1]<<" "<<10*backbone[k][2]<<"\n";
      }*/
    
  }

  double rx, ry, rz=0;
  rx=backbone[14][0]-backbone[2][0];
  ry=backbone[14][1]-backbone[2][1];
  rz=backbone[14][2]-backbone[2][2];
  double r2 = sqrt(rx*rx+ry*ry+rz*rz);

  print(counter, allpositions, r2, allenergies);
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
  beta[0]=0.400914;
  beta[1]=0.359021;
  beta[2]=0.32158;
  beta[3]=0.288422;
  beta[4]=0.25809;
  beta[5]=0.23084;
  beta[6]=0.20703;
}

void print(int counter, vector<double> & allpositions, double r2, vector<double> & allenergies){
  ofstream colvarfile;
  colvarfile.open("colvar.data", ios_base::app);

  for(int i=0; i<N_cvs; i++){
    colvarfile << allpositions[i];
    colvarfile << " ";
  }
  colvarfile << allenergies[0];
  colvarfile << " ";
  colvarfile << allenergies[1];
  colvarfile << " ";
  colvarfile << allenergies[2];
  colvarfile << " ";
  colvarfile << allenergies[3];
  colvarfile << " ";
  colvarfile << r2;
  colvarfile << "\n";

  colvarfile.close();

}
