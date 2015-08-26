#include <Two_d_grid.h> //includes the 2d grid class for handling everything you want to do with a 2d grid
#include <random>
#include <chrono>
#include <armadillo>

//g++ -O2 mc_dihedrals_pt.cpp -I . -I ~/code/armadillo-4.450.2/include -DARMA_USE_BLAS -DARMA_USE_LAPACK -DARMA_DONT_USE_WRAPPER -std=c++0x -lblas -llapack -o mc_dihedrals_pt

using namespace arma;
using namespace std;
int Nbiases=0;
int gridx=200;
int gridy=200;
int N_cvs=0;
double step_size=0.1;
vector <double> beta(8);
vector <double> crystal_dihedrals(50);
int Nbackbone=0;
int Nsteps=8000;
int Nsweeps=5000000;
vector <Two_d_grid> bias_grid;

void initialize(string);
void run_mc();
void print(int, vector<double> &, double, double, double, double,  vector<double> &);
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
  int Nreplicas=3;
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
  double accepted=0;
  double attempts=0;
  positions.resize(Nreplicas);
  
  for(int k=0; k<positions.size(); k++){
    positions[k].resize(N_cvs);
  }
  
  for(int i=0; i<=Nsweeps; i++){
    //cout << i << "\n"; 
    for(int k=0; k<Nreplicas; k++){
      //cout << k << "\n"; 
      for(int j=0; j<=Nsteps; j++){
	//cout << j << "\n"; 
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
    
    if(i%2==0){

      backbone(i,positions[0],energy);

    }
    swap_pair=distributionB(generator);
    double bf = (beta[swap_pair]-beta[swap_pair+1])*(energy[swap_pair+1]-energy[swap_pair]);
    
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
    }
  }
  
}

void backbone(int counter, vector<double> & allpositions, vector<double> & allenergies){
  vector< vector <double> > backbone;
  vector< double > bond_distances(Nbackbone+6);
  vector< double > bond_angles(Nbackbone+6);
  vector< double > dihedral_angles(Nbackbone+6);
  ofstream backbonefile;
  double dihedral_rmsd1=0;
  double dihedral_rmsd2=0;
  string bbstring ("backbone_atoms.xyz");
  backbonefile.open(bbstring,ios::out);
  backbone.resize(Nbackbone);
  for(int k=0; k<Nbackbone; k++){
    backbone[k].resize(3);
  }
  backbone[0][0]=0;backbone[0][1]=0;backbone[0][2]=0;
  bond_distances[0]=0.133386;
  for(int k=0; k<=N_cvs/2+1; k++){
    bond_distances[3*k]=0.147101;
    bond_distances[1+3*k]=0.1546;
    bond_distances[2+3*k]=0.133649;
  }
  bond_angles[0]=2.1957;
  for(int k=0; k<=N_cvs/2+1; k++){
    bond_angles[3*k]=2.0031;
    bond_angles[1+3*k]=2.0653;
    bond_angles[2+3*k]=2.17011;
  }
  int j1=0;

  for(int k=0; k<Nbackbone/3; k++){
    dihedral_angles[3*k+0]=allpositions[j1];j1++;
    dihedral_angles[3*k+1]=3.14159265;
    dihedral_angles[3*k+2]=allpositions[j1];j1++;
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
    
  }


  double rx, ry, rz=0;
  rx=backbone[28][0]-backbone[1][0];
  ry=backbone[28][1]-backbone[1][1];
  rz=backbone[28][2]-backbone[1][2];
  double r2 = sqrt(rx*rx+ry*ry+rz*rz);
  double meanx, meany, meanz=0;
  meanx=(backbone[1][0]+backbone[4][0]+backbone[7][0]+backbone[10][0]+backbone[13][0]+backbone[16][0]+backbone[19][0]+backbone[22][0]+backbone[25][0]+backbone[28][0])/10;
  meany=(backbone[1][1]+backbone[4][1]+backbone[7][1]+backbone[10][1]+backbone[13][1]+backbone[16][1]+backbone[19][1]+backbone[22][1]+backbone[25][1]+backbone[28][1])/10;
  meanz=(backbone[1][2]+backbone[4][2]+backbone[7][2]+backbone[10][2]+backbone[13][2]+backbone[16][2]+backbone[19][2]+backbone[22][2]+backbone[25][2]+backbone[28][2])/10;
  double rg = sqrt(pow(backbone[1][0]-meanx,2)+pow(backbone[1][1]-meany,2)+pow(backbone[1][2]-meanz,2));
  rg += sqrt(pow(backbone[4][0]-meanx,2)+pow(backbone[4][1]-meany,2)+pow(backbone[4][2]-meanz,2));
  rg += sqrt(pow(backbone[7][0]-meanx,2)+pow(backbone[7][1]-meany,2)+pow(backbone[7][2]-meanz,2));
  rg += sqrt(pow(backbone[10][0]-meanx,2)+pow(backbone[10][1]-meany,2)+pow(backbone[10][2]-meanz,2));
  rg += sqrt(pow(backbone[13][0]-meanx,2)+pow(backbone[13][1]-meany,2)+pow(backbone[13][2]-meanz,2));
  rg += sqrt(pow(backbone[16][0]-meanx,2)+pow(backbone[16][1]-meany,2)+pow(backbone[16][2]-meanz,2));
  rg += sqrt(pow(backbone[19][0]-meanx,2)+pow(backbone[19][1]-meany,2)+pow(backbone[19][2]-meanz,2));
  rg += sqrt(pow(backbone[22][0]-meanx,2)+pow(backbone[22][1]-meany,2)+pow(backbone[22][2]-meanz,2));
  rg += sqrt(pow(backbone[25][0]-meanx,2)+pow(backbone[25][1]-meany,2)+pow(backbone[25][2]-meanz,2));
  rg += sqrt(pow(backbone[28][0]-meanx,2)+pow(backbone[28][1]-meany,2)+pow(backbone[28][2]-meanz,2));
  rg = rg/10;
  
  rx=0;
  dihedral_rmsd2=0; dihedral_rmsd1=0;
  for( int k=0; k<N_cvs; k++){
    rx=allpositions[k]-crystal_dihedrals[k];
    if(rx>3.14159){
      rx=-3.14159+(rx-3.14159);
    }
    if(rx<-3.14159){
      rx=3.14159-(-3.14159-rx);
      }
    if(k<8){
      //if(k!=2){
	dihedral_rmsd1+=rx*rx;
	//}
      //dihedral_rmsd1+=0.5+0.5*cos(rx);
    }
    else{
      dihedral_rmsd2+=rx*rx;
      //dihedral_rmsd2+=0.5+0.5*cos(rx);
    }
  }
  dihedral_rmsd1=sqrt(dihedral_rmsd1/8);
  dihedral_rmsd2=sqrt(dihedral_rmsd2/10);
  //dihedral_rmsd1=dihedral_rmsd1;
  //dihedral_rmsd2=dihedral_rmsd2;
  print(counter, allpositions, r2, rg, dihedral_rmsd1, dihedral_rmsd2, allenergies);

  //cout << "hello 7\n";
  //cout << r2 <<" "<< bond_distances[3]<<" r2\n";
}


void initialize(string biasfname){
  Two_d_grid gridtemplate(gridx, gridy, 2*3.141592654, 2*3.141592654, -3.141592654, -3.141592654);
//Nbiases is the number of 2d_grids you want to make, gridx, and gridy are the number of gridpoints you are using
  for(int i=0; i<Nbiases; i++){
    bias_grid.push_back(gridtemplate);
    bias_grid[i].readfiles(i,biasfname);
    bias_grid[i].setindices();
  }

  Nbackbone=3*N_cvs/2 +2;
  beta[0]=0.353742*340/340;
  beta[1]=0.31684*340/340;
  beta[2]=0.28379*340/340;
  beta[3]=0.25419*340/340;
  beta[4]=0.22767*340/340;
  beta[5]=0.20392*340/340;
  beta[6]=0.18265*340/340;
  beta[7]=0.16359*340/340;
  ofstream colvarfile;
  //colvarfile.open("colvar_pt5rep.data", ios_base::app);
  colvarfile.open("colvar.data", ios_base::app);
  colvarfile <<"#! FIELDS time psi-1 phi-2 psi-2 phi-3 psi-3 phi-4 psi-4 phi-5 psi-5 phi-6 psi-6 phi-7 psi-7 phi-8 psi-8 phi-9 psi-9 phi-10 d1 rg dihedral_rmsd1 dihedral_rmsd2"<<"\n";
  colvarfile.close();
  crystal_dihedrals[0]=2.14;
  crystal_dihedrals[1]=-2.52;
  crystal_dihedrals[2]=1.98;
  crystal_dihedrals[3]=-1.58;
  crystal_dihedrals[4]=1.98;
  crystal_dihedrals[5]=-1.34;
  crystal_dihedrals[6]=-0.26;
  crystal_dihedrals[7]=-0.89011;
  crystal_dihedrals[8]=-1.1519;
  crystal_dihedrals[9]=-1.5882;
  crystal_dihedrals[10]=-0.36651;
  crystal_dihedrals[11]=1.4486;
  crystal_dihedrals[12]=0.078;
  crystal_dihedrals[13]=-1.186;
  crystal_dihedrals[14]=1.867;
  crystal_dihedrals[15]=-1.308;
  crystal_dihedrals[16]=2.356;
  crystal_dihedrals[17]=-2.635;

}

void print(int counter, vector<double> & allpositions, double r2, double rg, double dihedral_rmsd1, double dihedral_rmsd2, vector<double> & allenergies){
  ofstream colvarfile;
  //colvarfile.open("colvar_pt5rep.data", ios_base::app);
  colvarfile.open("colvar.data", ios_base::app);
  colvarfile << counter<<" ";
  for(int i=0; i<N_cvs; i++){
    colvarfile << allpositions[i];
    colvarfile << " ";
  }
  colvarfile << r2 <<" ";
  colvarfile << rg <<" ";
  colvarfile << dihedral_rmsd1 <<" ";
  colvarfile << dihedral_rmsd2 <<" ";
  colvarfile << "\n";

  colvarfile.close();

}
