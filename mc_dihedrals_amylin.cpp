#include <Two_d_grid.h> //includes the 2d grid class for handling everything you want to do with a 2d grid
#include <One_d_grid.h>
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
double pos_step_size=0.2;
double orientation_step_size=0.3;
vector <double> beta(8);
vector <double> crystal_dihedrals(50);
vector <double> alphahelix_dihedrals(50);
vector <double> betasheet_dihedrals(50);
vector <vector <double> > alphabeta;
vector <double> initial_dihedrals(50);
vector <vector <vector <vector <double> > > > atomistic_configurations;
int Nbackbone=0;
int Nreplicas=0;
int Nsteps=8000;
int Nsweeps=5000000;
//int Nsteps=1000;
//int Nsweeps=1;
int N_proteins=0;
double kappa=0;
double bessel=0;
double fold_weight=0;
//int Nsteps=0;
//int Nsweeps=0;

vector <Two_d_grid> bias_grid;
One_d_grid alphabetabias_helix(200,13,0);
One_d_grid alphabetabias_sheet(200,13,0);
Two_d_grid alphabeta_helix_distance(200,200,13,3.4,0,0.2);
Two_d_grid alphabeta_sheet_distance(200,200,13,3.4,0,0.2);
Two_d_grid angle_distance(200,200,3.141592654,3.4,0,0.2);

void initialize(string, string, string, string, string, string);
void run_mc();
void print(int, vector< vector <double> > &,  vector< vector <double> > &, vector <double> &, vector <double> &);
void backbone(int, vector <double> &, vector <double> &, vector <double> &, int, int);
int check_for_overlaps(int);
void print_trajectory(int);

string outfile;
string trajectoryfile;
string xyzfile;

int main(int argc, char *argv[]){
  string biasfilename;
  string ab_helix_filename;
  string ab_sheet_filename;
  string ab_helix_distance_filename;
  string ab_sheet_distance_filename;
  string angle_distance_filename;
  
  if(argc>1){
    N_cvs = atoi( argv[1]);
    Nbiases = atoi (argv[2]);
    N_proteins = atoi (argv[3]);
    biasfilename = argv[4];
    ab_helix_filename = argv[5];
    ab_sheet_filename = argv[6];
    ab_helix_distance_filename = argv[7];
    ab_sheet_distance_filename = argv[8];
    angle_distance_filename = argv[9];
    outfile=argv[10];
    xyzfile=argv[11];
    kappa = atof(argv[12]);
    bessel = atof(argv[13]);
  }
  
  initialize(biasfilename,ab_helix_filename,ab_sheet_filename,ab_helix_distance_filename,ab_sheet_distance_filename,angle_distance_filename);
  
  
  run_mc();
  
  return 0; 
}

void run_mc(){
  vector <vector <vector <double> > > dihedral_angles(Nreplicas, vector <vector <double> >(N_proteins, vector <double>(N_cvs))); //index order: 1st corresponds to replica, 2nd to protein, 3rd to dihedral 
  vector <vector <vector <double> > > positions(Nreplicas, vector <vector <double> >(N_proteins, vector <double>(3)));
  vector <vector <vector <double> > > orientations(Nreplicas, vector <vector <double> >(N_proteins, vector <double>(3)));
  vector <vector <double> > ab_alphahelix(Nreplicas, vector <double>(N_proteins));
  vector <vector <double> > ab_betasheet(Nreplicas, vector <double>(N_proteins));
  vector <vector <double> > pofs_helix(Nreplicas, vector <double>(N_proteins));
  vector <vector <double> > pofs_sheet(Nreplicas, vector <double>(N_proteins));
  vector<double> pos_temp(3);
  vector<double> orientation_temp(3);
  vector<double> pos_old(3);
  vector<double> orientation_old(3);
  double angle_old=0;
  double angle_temp=0;
  double energy_delta=0;
  double alphabeta1new=0;
  double ab_alphahelix_new=0;
  double ab_betasheet_new=0;
  double alphabeta2new=0;
  vector <double> energy(Nreplicas,0);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  uniform_real_distribution<double> real_distribution(0.0,1.0);
  uniform_int_distribution<int> distributionA(0,N_cvs-1);
  uniform_int_distribution<int> distributionB(0,Nreplicas-2);
  uniform_int_distribution<int> distributionC(0,N_proteins-1);
  int dihedral_label=0;
  int protein_label=0;
  int swap_pair=0;
  double xmax=3.141592654;
  double xmin=-3.14159264;
  double accepted=0;
  double pofsold=0;
  double pofsnew_ah=0;
  double pofsnew_bs=0;
  double pofs0=0;
  double attempts=0;
  double increment=0;
  double rx=0;
  double rx2=0;
  double pofsinc=0;
  double r1=0;
  double r2=0;
  double r3=0;
  double r0=0;
  double rsq=0;
  double d1=0;
  double d2=0;
  double d3=0;
  double d0=0;
  double flat_norm=0;
  int outofbounds=0;
  int Noverlaps=0;
  //flat_norm = 1/(pow(2*3.1415926,2*N_cvs)*3.4*3.1415926);
  flat_norm = 1/(3.4*3.1415926);

  /*for(int k=0; k<positions.size(); k++){
    for(int j=0; j<N_cvs; j++){
      if(j<15){
	positions[k][j]=crystal_dihedrals[j]+0.3;
	positions[k][j]=initial_dihedrals[j];
      }
      else{
	positions[k][j]=crystal_dihedrals[j];
	positions[k][j]=initial_dihedrals[j];
      }
    }
    }*/
  
  
  for(int k=0; k<positions.size(); k++){
    for(int j=0; j<positions[k].size(); j++){
      
      outofbounds=1;
      while(outofbounds !=0){
	positions[k][j][0]=2.07*real_distribution(generator);
	positions[k][j][1]=2.07*real_distribution(generator);
	positions[k][j][2]=2.07*real_distribution(generator);
	outofbounds=0;
	for(int i1=0; i1<j; i1++){
	  r2 = sqrt( pow(positions[k][j][0]-positions[k][i1][0],2)+pow(positions[k][j][1]-positions[k][i1][1],2)+ pow(positions[k][j][2]-positions[k][i1][2],2) );
	  //cout << r2 << "\n";
	  if(r2<0.2 || r2>3.6){
	    outofbounds+=1;
	  }
	}
      }
      
      //cout << r2 << " a\n";
      //cout << positions[k][j][0] << " " << positions[k][j][1] << " " << positions[k][j][2] << "\n";
      
      orientations[k][j][0]=real_distribution(generator);
      orientations[k][j][1]=real_distribution(generator);
      orientations[k][j][2]=real_distribution(generator);
      double mag= sqrt(pow(orientations[k][j][0],2) + pow(orientations[k][j][1],2) +pow(orientations[k][j][2],2));
      orientations[k][j][0] /= mag;
      orientations[k][j][1] /= mag;
      orientations[k][j][2] /= mag;
      ab_alphahelix[k][j]=0;
      ab_betasheet[k][j]=0;
      for(int j1=0; j1<dihedral_angles[k][j].size(); j1++){
	//dihedral_angles[k][j][j1]=-3.14+real_distribution(generator)*2*3.14;
	//dihedral_angles[k][j][j1]=alphahelix_dihedrals[j1]+3.14/2;
	dihedral_angles[k][j][j1]=alphahelix_dihedrals[j1];
	double angle_delta = dihedral_angles[k][j][j1]-alphahelix_dihedrals[j1];
	ab_alphahelix[k][j] += 0.5+0.5*cos(angle_delta);
	angle_delta = dihedral_angles[k][j][j1]-betasheet_dihedrals[j1];
	ab_betasheet[k][j] += 0.5+0.5*cos(angle_delta);
      }
      //cout << ab_alphahelix[k][j] << " ab ah \n";
      backbone(0, dihedral_angles[k][j], orientations[k][j], positions[k][j], k, j);
      
      //cout << ab_alphahelix[k][j] << " ab ah \n";
    }
  
  }
  print_trajectory(0);

  for(int k=0; k<pofs_helix.size(); k++){
    for(int j=0; j<pofs_helix.size(); j++){
      pofs_helix[k][j] = 1;
      pofs_sheet[k][j] = 1;
      for(int i=0; i<N_cvs; i++){
	pofs_helix[k][j] *= exp(kappa*cos(dihedral_angles[k][j][i]-alphahelix_dihedrals[i]))/(bessel);
	pofs_sheet[k][j] *= exp(kappa*cos(dihedral_angles[k][j][i]-betasheet_dihedrals[i]))/(bessel);
	
      }
      
    }

  }
  

  
  for(int i=0; i<=Nsweeps; i++){
    //cout << i << "\n"; 
    for(int k=0; k<Nreplicas; k++){//Nreplicas don't forget to change back
      //cout << k << "\n"; 
      for(int j=0; j<=Nsteps; j++){
	//cout << j << "\n"; 
	dihedral_label = distributionA(generator);
	protein_label = distributionC(generator);
	angle_temp = dihedral_angles[k][protein_label][dihedral_label] + (real_distribution(generator)-0.5)*step_size;
	pos_temp[0] = positions[k][protein_label][0]+(real_distribution(generator)-0.5)*pos_step_size;
	pos_temp[1] = positions[k][protein_label][1]+(real_distribution(generator)-0.5)*pos_step_size;
	pos_temp[2] = positions[k][protein_label][2]+(real_distribution(generator)-0.5)*pos_step_size;

	pos_old[0] = positions[k][protein_label][0];
	pos_old[1] = positions[k][protein_label][1];
	pos_old[2] = positions[k][protein_label][2];
	
	if(real_distribution(generator)<0.2){
	  angle_temp=angle_temp+3.14;
	}
	
	if(angle_temp>xmax-0.01){
	  angle_temp=xmin+(angle_temp-xmax);
	}
	else if(angle_temp<xmin){
	  angle_temp=xmax-(xmin-angle_temp)-0.001;
	}
	rsq=2;
	while(rsq >= 1){
	  r1= 1-2*real_distribution(generator);
	  r2= 1-2*real_distribution(generator);   
	  rsq=r1*r1 + r2*r2;
	}
	r0=2*sqrt(1-rsq);
	r1=orientation_step_size*r1*r0;
	r2=orientation_step_size*r2*r0;
	r3=orientation_step_size*(1-2*rsq);
	d1=orientations[k][protein_label][0]+r1;
	d2=orientations[k][protein_label][1]+r2;
	d3=orientations[k][protein_label][2]+r3;
	d0=d1*d1 + d2*d2 + d3*d3;
	d0=sqrt(d0);
	orientation_temp[0]=d1/d0;
	orientation_temp[1]=d2/d0;
	orientation_temp[2]=d3/d0;
	
	orientation_old[0]=orientations[k][protein_label][0];
	orientation_old[1]=orientations[k][protein_label][1];
	orientation_old[2]=orientations[k][protein_label][2];
	
	angle_old=dihedral_angles[k][protein_label][dihedral_label];
	//pofsnew_ah = pofs_helix[k][j]*exp(kappa*cos(angle_temp-alphahelix_dihedrals[dihedral_label]))/(bessel);
	//pofsnew_ah /= exp(kappa*cos(dihedral_angles[k][protein_label][dihedral_label]-alphahelix_dihedrals[dihedral_label]))/(bessel);
	//pofsnew_bs = pofs_sheet[k][j]*exp(kappa*cos(angle_temp-betasheet_dihedrals[dihedral_label]))/(bessel);
	//pofsnew_bs /= exp(kappa*cos(dihedral_angles[k][protein_label][dihedral_label]-betasheet_dihedrals[dihedral_label]))/(bessel);
	//pofsnew_ah=1;
	//pofsnew_bs=1;
	//cout << pofsinc << " " << pofsold << " " << pofsnew << " \n";
	
	ab_alphahelix_new=0;
	ab_betasheet_new=0;
	for(int j1=0; j1<N_cvs; j1++){
	  if(j1==dihedral_label){
	    rx=angle_temp-alphahelix_dihedrals[dihedral_label];
	    rx2=angle_temp-betasheet_dihedrals[dihedral_label];
	  }
	  else{
	    rx=dihedral_angles[k][protein_label][j1]-alphahelix_dihedrals[j1];
	    rx2=dihedral_angles[k][protein_label][j1]-betasheet_dihedrals[j1];
	  }
	  ab_alphahelix_new+=0.5+0.5*cos(rx);
	  ab_betasheet_new+=0.5+0.5*cos(rx2);
	}
	pofsnew_ah = exp(kappa*(2*ab_alphahelix_new-N_cvs))/(pow(bessel,N_cvs));
	pofsnew_bs = exp(kappa*(2*ab_betasheet_new-N_cvs))/(pow(bessel,N_cvs));	
	
	energy_delta=0;
	
	if(dihedral_label>0 && dihedral_label< N_cvs-1){
	  energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(angle_temp, dihedral_angles[k][protein_label][dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label],dihedral_angles[k][protein_label][dihedral_label+1]);
	  energy_delta += bias_grid[dihedral_label-1].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label-1], angle_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label-1],dihedral_angles[k][protein_label][dihedral_label]);

	}
	
	else if(dihedral_label==0){	
	  energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(angle_temp, dihedral_angles[k][protein_label][dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label],dihedral_angles[k][protein_label][dihedral_label+1]);

	}
	else{
	  energy_delta = bias_grid[dihedral_label-1].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label-1], angle_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label-1],dihedral_angles[k][protein_label][dihedral_label]);

	}
	//cout << energy_delta << " a\n";
	energy_delta += alphabetabias_helix.getvalue_linearinterpolation(ab_alphahelix_new)-alphabetabias_helix.getvalue_linearinterpolation(ab_alphahelix[k][protein_label]);
	//cout << energy_delta << " b\n";
	energy_delta += alphabetabias_sheet.getvalue_linearinterpolation(ab_betasheet_new)-alphabetabias_sheet.getvalue_linearinterpolation(ab_betasheet[k][protein_label]);
	//cout << energy_delta << " c\n";
	//cout << "a1 "<< energy_delta << " " << angle_temp << " " << dihedral_angles[k][protein_label][dihedral_label] << " " << ab_alphahelix_new << " " << ab_betasheet_new << "\n";
	outofbounds=0;
	
	for(int i1=0; i1< N_proteins; i1++){
	  if(i1 != protein_label){
	    double ip_r=0;
	    double ip_r_new=0;
	    double interproteinx = positions[k][protein_label][0]-positions[k][i1][0];
	    double interproteiny = positions[k][protein_label][1]-positions[k][i1][1];
	    double interproteinz = positions[k][protein_label][2]-positions[k][i1][2];
	    double interproteinx_new = pos_temp[0]-positions[k][i1][0];
	    double interproteiny_new = pos_temp[1]-positions[k][i1][1];
	    double interproteinz_new = pos_temp[2]-positions[k][i1][2];
	    double pofs_new=0;
	    ip_r = pow(interproteinx*interproteinx + interproteiny*interproteiny + interproteinz*interproteinz, 0.5);
	    ip_r_new = pow(interproteinx_new*interproteinx_new + interproteiny_new*interproteiny_new + interproteinz_new*interproteinz_new, 0.5);
	    //cout << "a2" << " " << ip_r << " " << ip_r_new << " " << ab_alphahelix_new << " " << ab_alphahelix[k][protein_label] << "\n";

	    energy_delta += alphabeta_helix_distance.getvalue_linearinterpolation(ab_alphahelix_new,ip_r_new)-alphabeta_helix_distance.getvalue_linearinterpolation(ab_alphahelix[k][protein_label],ip_r);
	    energy_delta += alphabeta_helix_distance.getvalue_linearinterpolation(ab_alphahelix[k][i1],ip_r_new)-alphabeta_helix_distance.getvalue_linearinterpolation(ab_alphahelix[k][i1],ip_r);
	    //cout << "loop 1 " <<  energy_delta << "\n";	    
	    energy_delta += alphabeta_sheet_distance.getvalue_linearinterpolation(ab_betasheet_new,ip_r_new)-alphabeta_sheet_distance.getvalue_linearinterpolation(ab_betasheet[k][protein_label],ip_r);
	    energy_delta += alphabeta_sheet_distance.getvalue_linearinterpolation(ab_betasheet[k][i1],ip_r_new)-alphabeta_sheet_distance.getvalue_linearinterpolation(ab_betasheet[k][i1],ip_r);
	    //cout << "loop 2 "<<  energy_delta << "\n";
	    double angle_new = acos(orientation_temp[0]*orientations[k][i1][0] + orientation_temp[1]*orientations[k][i1][1] + orientation_temp[2]*orientations[k][i1][2]);
	    double angle = acos(orientations[k][protein_label][0]*orientations[k][i1][0] + orientations[k][protein_label][1]*orientations[k][i1][1] + orientations[k][protein_label][2]*orientations[k][i1][2]);
	    
	    //cout << "b" << " " << angle << " " << angle_new << " " << orientations[k][protein_label][0]<< "\n";
	    //angle=angle-0.1;
	    energy_delta += angle_distance.getvalue_linearinterpolation(angle_new,ip_r_new)-angle_distance.getvalue_linearinterpolation(angle,ip_r);
	    //cout << "b1" << energy_delta << "\n";
	    double pofsnew = pofsnew_ah*pofs_helix[k][i1]*exp(-(ip_r_new-0.5)*(ip_r_new-0.5)/(2*0.3*0.3))*exp(2.5*cos(angle_new-3.14159))/(0.632682*10.3353*ip_r_new*ip_r_new*sin(angle_new)); 
	    //double pofsnew = exp(-(ip_r_new-0.5)*(ip_r_new-0.5)/(2*0.3*0.3))*exp(2.5*cos(angle_new-3.14159))/(0.632682*10.3353*ip_r_new*ip_r_new*sin(angle_new));
	    //double pofsnew = pofsnew_ah*pofs_helix[k][i1]; 
	    
	    pofsnew += pofsnew_bs*pofs_sheet[k][i1]*exp(-(ip_r_new-0.5)*(ip_r_new-0.5)/(2*0.3*0.3))*exp(2.5*cos(angle_new-3.14159))/(0.632682*10.3353*ip_r_new*ip_r_new*sin(angle_new));
	    //pofsnew += exp(-(ip_r_new-0.5)*(ip_r_new-0.5)/(2*0.3*0.3))*exp(2.5*cos(angle_new-3.14159))/(0.632682*10.3353*ip_r_new*ip_r_new*sin(angle_new));	    
	    //pofsnew += pofsnew_bs*pofs_sheet[k][i1];
	    
	    pofsnew += 1/(ip_r_new*ip_r_new*sin(angle_new)*3.1415926*3.4);
	    //pofsnew += 1;
	    
	    double pofsold = pofs_helix[k][protein_label]*pofs_helix[k][i1]*exp(-(ip_r-0.5)*(ip_r-0.5)/(2*0.3*0.3))*exp(2.5*cos(angle-3.14159))/(0.632682*10.3353*ip_r*ip_r*sin(angle));
	    pofsold += pofs_sheet[k][protein_label]*pofs_sheet[k][i1]*exp(-(ip_r-0.5)*(ip_r-0.5)/(2*0.3*0.3))*exp(2.5*cos(angle-3.14159))/(0.632682*10.3353*ip_r*ip_r*sin(angle));
	    //double pofsold = exp(-(ip_r-0.5)*(ip_r-0.5)/(2*0.3*0.3))*exp(2.5*cos(angle-3.14159))/(0.632682*10.3353*ip_r*ip_r*sin(angle));
	    //pofsold += exp(-(ip_r-0.5)*(ip_r-0.5)/(2*0.3*0.3))*exp(2.5*cos(angle-3.14159))/(0.632682*10.3353*ip_r*ip_r*sin(angle));
	    
	    //double pofsold = pofs_helix[k][protein_label]*pofs_helix[k][i1];
	    //pofsold += pofs_sheet[k][protein_label]*pofs_sheet[k][i1];
	    
	    pofsold += 1/(ip_r*ip_r*sin(angle)*3.1415926*3.4);
	    //pofsold += 1;
	    //cout << pofsold << " " << pofsnew<< " "<< pofsnew_ah<<" "<<pofsnew_bs<< " " << pofs_helix[k][protein_label]<< " "<< pofs_sheet[k][protein_label]<< " c\n";
	    //cout << energy_delta << " a\n";
	    //energy_delta=0;
	    energy_delta += (log(pofsold)-log(pofsnew))/beta[0];
	    //cout << energy_delta << " b\n";
	  
	    if(ip_r_new>3.6 || ip_r_new<0.2){
	      outofbounds+=1;
	      //cout << outofbounds << " outofbounds \n";
	    }
	    if(angle_new < 1.3 && angle_new>0.2 && ip_r_new<0.8){
	      outofbounds+=1;
	    }
	    //cout << energy_delta << " " << outofbounds << "\n";
	  }

	  
	}
	//cout << energy_delta << " d\n";
	
	
	if((energy_delta<0 && outofbounds==0) || (real_distribution(generator)<=exp(-beta[k]*energy_delta) && outofbounds==0)){
	  energy[k] += energy_delta;
	  dihedral_angles[k][protein_label][dihedral_label]=angle_temp;
	  positions[k][protein_label][0]=pos_temp[0];
	  positions[k][protein_label][1]=pos_temp[1];
	  positions[k][protein_label][2]=pos_temp[2];
	  //cout << outofbounds << " accepted\n";
	  orientations[k][protein_label][0]=orientation_temp[0];
	  orientations[k][protein_label][1]=orientation_temp[1];
	  orientations[k][protein_label][2]=orientation_temp[2];	  
	  backbone(0, dihedral_angles[k][protein_label], orientations[k][protein_label], positions[k][protein_label], k, protein_label);	  
	  
	  Noverlaps = check_for_overlaps(0);
	  if(Noverlaps !=0){
	    energy[k] -= energy_delta;
	    dihedral_angles[k][protein_label][dihedral_label]=angle_old;
	    positions[k][protein_label][0]=pos_old[0];
	    positions[k][protein_label][1]=pos_old[1];
	    positions[k][protein_label][2]=pos_old[2];
	    //cout << outofbounds << " accepted\n";
	    orientations[k][protein_label][0]=orientation_old[0];
	    orientations[k][protein_label][1]=orientation_old[1];
	    orientations[k][protein_label][2]=orientation_old[2];	  
	    backbone(0, dihedral_angles[k][protein_label], orientations[k][protein_label], positions[k][protein_label], k, protein_label);
	  }
	  if(Noverlaps==0){
	    ab_alphahelix[k][protein_label]=ab_alphahelix_new;
	    ab_betasheet[k][protein_label]=ab_betasheet_new;
	  
	    pofs_helix[k][protein_label]=pofsnew_ah;
	    pofs_sheet[k][protein_label]=pofsnew_bs;
	  }
	  //cout << pofs_helix[k][protein_label] << " helix\n";
	  //cout << pofs_sheet[k][protein_label] << " sheet\n";
	}
	else{
	  //cout << "rejected\n";
	  angle_temp=angle_temp;
	  if(k==0){
	    pofs0=pofsold;
	  }
	}
      }
    }
    
    if(i%2==0){
      for(int j=0; j<N_proteins; j++){
	backbone(i, dihedral_angles[0][j], orientations[0][j], positions[0][j], 0, j);
      }
      print(i, positions[0], orientations[0], ab_alphahelix[0], ab_betasheet[0]);
      print_trajectory(i);
    }
    
    //swap_pair=distributionB(generator);
    //double bf = (beta[swap_pair]-beta[swap_pair+1])*(energy[swap_pair+1]-energy[swap_pair]);
    
    /*if(bf>0 || real_distribution(generator)<exp(bf)){
      for(int k=0; k<N_cvs; k++){
	double tmp1=dihedral_angles[swap_pair][k];
	double tmp2=dihedral_angles[swap_pair+1][k];
	dihedral_angles[swap_pair+1][k]=tmp1;
	dihedral_angles[swap_pair][k]=tmp2;	
      }
      for(int k=0; k<2; k++){
	double tmp1=alphabeta[swap_pair][k];
	double tmp2=alphabeta[swap_pair+1][k];
	alphabeta[swap_pair+1][k]=tmp1;
	alphabeta[swap_pair][k]=tmp2;		
      }
      double tmp1=energy[swap_pair];
      double tmp2=energy[swap_pair+1];
      energy[swap_pair+1]=tmp1;
      energy[swap_pair]=tmp2;
      }*/
  }

}


void print_trajectory(int timestep){
  ofstream backbonefile;

  backbonefile.open(xyzfile,ios::app);
  
  int Noverlaps = check_for_overlaps(0);
  //cout << Noverlaps << " no\n";
  if(Noverlaps==0){
    //cout << Noverlaps << " no2 printing\n";
    backbonefile <<N_proteins*Nbackbone << "\n\n";
    for(int i=0; i<N_proteins; i++){
      for(int j=0; j<Nbackbone; j++){
	if( (j-1) % 3 == 0){
	  backbonefile << "N " << 10*(atomistic_configurations[0][i][j][0]-atomistic_configurations[0][0][8][0]) << " " << 10*(atomistic_configurations[0][i][j][1]-atomistic_configurations[0][0][8][1]) << " " << 10*(atomistic_configurations[0][i][j][2]-atomistic_configurations[0][0][8][2]) << "\n";
	}
	else{
	  backbonefile << "C " << 10*(atomistic_configurations[0][i][j][0]-atomistic_configurations[0][0][8][0]) << " " << 10*(atomistic_configurations[0][i][j][1]-atomistic_configurations[0][0][8][1]) << " " << 10*(atomistic_configurations[0][i][j][2]-atomistic_configurations[0][0][8][2]) << "\n";
	}
      }
    }
  }
  //backbonefile.close();
}


int check_for_overlaps(int protein_number){
  int Noverlaps=0;
  for(int i=0; i<N_proteins; i++){
    for(int i1=i+1; i1<N_proteins; i1++){
      for(int j=0; j<Nbackbone; j++){
	for(int j1=0; j1<Nbackbone; j1++){
	  double rx = atomistic_configurations[0][i][j][0]-atomistic_configurations[0][i1][j1][0];
	  double ry = atomistic_configurations[0][i][j][1]-atomistic_configurations[0][i1][j1][1];
	  double rz = atomistic_configurations[0][i][j][2]-atomistic_configurations[0][i1][j1][2];
	  double r2 = sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));
	  //cout << i << " "<< i1 << " "<< r2<< "\n";
	  if(r2<.2){
	    Noverlaps+=1;
	  }
	}
      }
    }
  }
  //cout << Noverlaps << "overlaps\n";
  return Noverlaps;
}


void backbone(int counter, vector<double> & allpositions, vector <double> & orientation, vector <double> & center_position, int replica_number, int protein_number){
  vector< double > bond_distances(Nbackbone+8);
  vector< double > bond_angles(Nbackbone+8);
  vector< double > dihedral_angles(Nbackbone+8);
  ofstream backbonefile;
  double dihedral_alphabeta1=0;
  double dihedral_alphabeta2=0;
  double alphahelix_ab=0;
  double betastrand_ab=0;
  string bbstring ("backbone_atoms.xyz");
  //cout << "hello\n";
  bond_distances[0]=0.133386;
  for(int k=0; k<=N_cvs/2+2; k++){
    bond_distances[3*k]=0.147101;
    bond_distances[1+3*k]=0.1546;
    bond_distances[2+3*k]=0.133649;
  }
  
  bond_angles[0]=2.1957;
  for(int k=0; k<=N_cvs/2+2; k++){
    bond_angles[3*k+1]=2.0031;
    bond_angles[2+3*k]=2.0653;
    bond_angles[3*k]=2.17011;
  }
  
  int j1=0;
  //cout << "hello\n";  
  //cout << allpositions.size() << " y\n";
  for(int k=0; k<Nbackbone/3-1; k++){
    //cout << k << " \n";
    dihedral_angles[3*k+0]=allpositions[j1];j1++;
    dihedral_angles[3*k+1]=allpositions[j1];j1++;
    dihedral_angles[3*k+2]= 3.14159265;
    //cout << j1 << " c\n";
  }
  
  //cout << atomistic_configurations[replica_number][protein_number][2][0] << " ok\n";  
  vec terminal_vector(4);
  mat transformation_matrix(4,4,fill::eye);
  mat new_matrix(4,4,fill::eye);
  new_matrix(0,0)=-1;
  new_matrix(2,2)=-1;
  new_matrix(0,3)=-bond_distances[0];
  terminal_vector(0)=0; terminal_vector(1)=0; terminal_vector(2)=0; terminal_vector(3)=1;
  
  atomistic_configurations[replica_number][protein_number][0][0]=0;atomistic_configurations[replica_number][protein_number][0][1]=0; atomistic_configurations[replica_number][protein_number][0][2]=0;
  
  atomistic_configurations[replica_number][protein_number][1][0]=-bond_distances[0];atomistic_configurations[replica_number][protein_number][1][1]=0; atomistic_configurations[replica_number][protein_number][1][2]=0;
  //cout << atomistic_configurations[replica_number][protein_number][2][0] << " \n";
  transformation_matrix=transformation_matrix*new_matrix;
  new_matrix(0,0) = -cos(bond_angles[0]); new_matrix(0,1) = -sin(bond_angles[0]); new_matrix(0,2)=0; new_matrix(0,3)= -bond_distances[1]*cos(bond_angles[0]);
  new_matrix(1,0) = sin(bond_angles[0]); new_matrix(1,1) = -cos(bond_angles[0]); new_matrix(1,2)=0; new_matrix(1,3)= bond_distances[1]*sin(bond_angles[0]);
  new_matrix(2,2)=1;
  //cout << atomistic_configurations[replica_number][protein_number][2][0] << " \n";  
  transformation_matrix=transformation_matrix*new_matrix;
  atomistic_configurations[replica_number][protein_number][2][0]=transformation_matrix(0,3);
  atomistic_configurations[replica_number][protein_number][2][1]=transformation_matrix(1,3);
  atomistic_configurations[replica_number][protein_number][2][2]=transformation_matrix(2,3);
  //cout << atomistic_configurations[replica_number][protein_number][2][0] << " \n";
  //backbonefile <<"N "<<10*backbone[0][0] <<" "<<10*backbone[0][1]<<" "<<10*backbone[0][2]<<"\n";
  //backbonefile <<"C "<<10*backbone[1][0] <<" "<<10*backbone[1][1]<<" "<<10*backbone[1][2]<<"\n";
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
    
    atomistic_configurations[replica_number][protein_number][k][0]=transformation_matrix(0,3);
    atomistic_configurations[replica_number][protein_number][k][1]=transformation_matrix(1,3);
    atomistic_configurations[replica_number][protein_number][k][2]=transformation_matrix(2,3);
    
    //if(k%3==0){
      //backbonefile <<"N " <<10*backbone[k][0] <<" "<<10*backbone[k][1]<<" "<<10*backbone[k][2]<<"\n";
    //}
    //else{
      //backbonefile <<"C " <<10*backbone[k][0] <<" "<<10*backbone[k][1]<<" "<<10*backbone[k][2]<<"\n";
      //}
  
  }
  vec orientation_vec(3);
  double refx=atomistic_configurations[replica_number][protein_number][2][0];   double refy=atomistic_configurations[replica_number][protein_number][2][1];   double refz=atomistic_configurations[replica_number][protein_number][2][2];
  
  for(int k=0; k<Nbackbone; k++){
    atomistic_configurations[replica_number][protein_number][k][0] -= refx;
    atomistic_configurations[replica_number][protein_number][k][1] -= refy;
    atomistic_configurations[replica_number][protein_number][k][2] -= refz;
  }
  orientation_vec(0) = atomistic_configurations[replica_number][protein_number][20][0]-atomistic_configurations[replica_number][protein_number][2][0];
  orientation_vec(1) = atomistic_configurations[replica_number][protein_number][20][1]-atomistic_configurations[replica_number][protein_number][2][1];  
  orientation_vec(2) = atomistic_configurations[replica_number][protein_number][20][2]-atomistic_configurations[replica_number][protein_number][2][2];
  refx= sqrt(pow(orientation_vec(0),2)+pow(orientation_vec(1),2)+pow(orientation_vec(2),2));
  orientation_vec(0) /= refx; orientation_vec(1) /= refx;   orientation_vec(2) /= refx;
  vec target_orientation(3);
  target_orientation(0)=orientation[0]; target_orientation(1)=orientation[1]; target_orientation(2)=orientation[2];
  vec cross_product(3);
  cross_product = cross(orientation_vec,target_orientation);
  refx = sqrt(pow(cross_product(0),2) + pow(cross_product(1),2) + pow(cross_product(2),2) );
  double theta = asin(refx);
  
  double sin_theta= sin(theta);   double cos_theta= cos(theta);
  
  mat rotation_matrix(3,3);
  orientation_vec(0) = cross_product(0)/refx;   orientation_vec(1) = cross_product(1)/refx;   orientation_vec(2) = cross_product(2)/refx;
  
  //cout << refx << " "<< orientation_vec(0) << " "<< orientation_vec(1) << " "<< orientation_vec(2)<< " " << sin_theta<< " " << cos_theta << " hey \n";
  
  //cout << sqrt(pow(orientation_vec(0),2) + pow(orientation_vec(1),2) + pow(orientation_vec(2),2) )<< " mag \n";
  rotation_matrix(0,0)= cos_theta + orientation_vec(0)*orientation_vec(0)*(1-cos_theta); rotation_matrix(0,1) = orientation_vec(0)*orientation_vec(1)*(1-cos_theta) - orientation_vec(2)*sin_theta; rotation_matrix(0,2)= orientation_vec(0)*orientation_vec(2)*(1-cos_theta) + orientation_vec(1)*sin_theta;
  rotation_matrix(1,0) = orientation_vec(0)*orientation_vec(1)*(1-cos_theta) + orientation_vec(2)*sin_theta; rotation_matrix(1,1)=cos_theta + orientation_vec(1)*orientation_vec(1)*(1-cos_theta); rotation_matrix(1,2)= orientation_vec(1)*orientation_vec(2)*(1-cos_theta) - orientation_vec(0)*sin_theta;
  rotation_matrix(2,0) = orientation_vec(0)*orientation_vec(2)*(1-cos_theta) - orientation_vec(1)*sin_theta; rotation_matrix(2,1) = orientation_vec(1)*orientation_vec(2)*(1-cos_theta) + orientation_vec(0)*sin_theta; rotation_matrix(2,2)= cos_theta + orientation_vec(2)*orientation_vec(2)*(1-cos_theta);
  
  vec unrotated_position(3);
  vec rotated_position(3);
  
  for(int k=0; k<Nbackbone; k++){
    unrotated_position(0)= atomistic_configurations[replica_number][protein_number][k][0]; unrotated_position(1)= atomistic_configurations[replica_number][protein_number][k][1]; unrotated_position(2)= atomistic_configurations[replica_number][protein_number][k][2];
    rotated_position = rotation_matrix*unrotated_position;
    //cout << sqrt(pow(rotated_position(0),2) + pow(rotated_position(1),2) + pow(rotated_position(2),2) )<< " " <<sqrt(pow(unrotated_position(0),2) + pow(unrotated_position(1),2) + pow(unrotated_position(2),2) )<< "\n";
    atomistic_configurations[replica_number][protein_number][k][0]=rotated_position(0); atomistic_configurations[replica_number][protein_number][k][1]=rotated_position(1); atomistic_configurations[replica_number][protein_number][k][2]=rotated_position(2);
    }
  refx = center_position[0] - atomistic_configurations[replica_number][protein_number][8][0]; refy = center_position[1] - atomistic_configurations[replica_number][protein_number][8][1]; refz = center_position[2] - atomistic_configurations[replica_number][protein_number][8][2];
  
  //cout << refx << " " << refy << " " << refz << "\n";
  //cout << atomistic_configurations[replica_number][protein_number][2][0] << " " << atomistic_configurations[replica_number][protein_number][2][1] << " " << atomistic_configurations[replica_number][protein_number][2][2] << "a\n";
  for(int k=0; k<Nbackbone; k++){
    atomistic_configurations[replica_number][protein_number][k][0] += refx; atomistic_configurations[replica_number][protein_number][k][1] += refy; atomistic_configurations[replica_number][protein_number][k][2] += refz;
  }
  //cout << atomistic_configurations[replica_number][protein_number][2][0] << " " << atomistic_configurations[replica_number][protein_number][2][1] << " " << atomistic_configurations[replica_number][protein_number][2][2] << "b\n";
  /*
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
  int noverlaps=0;
  for(int k=0; k<28; k++){
    for(int k1=0; k1<k-1; k1++){
      rx=sqrt( pow(backbone[k][0]-backbone[k1][0],2) +pow(backbone[k][1]-backbone[k1][1],2) +pow(backbone[k][2]-backbone[k1][2],2));
      //cout <<rx << "\n";
      if(rx < 0.17){
	noverlaps+=1;
	//cout << k << " " << k1 << "\n";
      }
    }
  }
  //cout << noverlaps << "\n";  
  
  rx=0;
  dihedral_alphabeta2=0; dihedral_alphabeta1=0;
  alphahelix_ab=0;
  betastrand_ab=0;
  for( int k=0; k<N_cvs; k++){
    rx=allpositions[k]-crystal_dihedrals[k];
    if(rx>3.14159){
      rx=-3.14159+(rx-3.14159);
    }
    if(rx<-3.14159){
      rx=3.14159-(-3.14159-rx);
    }
    if(k<8){
      //dihedral_alphabeta1+=rx*rx;
      dihedral_alphabeta1+=0.5+0.5*cos(rx);
    }
    else{
      //dihedral_alphabeta2+=rx*rx;
      dihedral_alphabeta2+=0.5+0.5*cos(rx);
    }
    rx=allpositions[k]-alphahelix_dihedrals[k];
    if(rx>3.14159){
      rx=-3.14159+(rx-3.14159);
    }
    if(rx<-3.14159){
      rx=3.14159-(-3.14159-rx);
    }
    alphahelix_ab+=0.5+0.5*cos(rx);

    rx=allpositions[k]-betastrand_dihedrals[k];
    if(rx>3.14159){
      rx=-3.14159+(rx-3.14159);
    }
    if(rx<-3.14159){
      rx=3.14159-(-3.14159-rx);
    }
    betastrand_ab+=0.5+0.5*cos(rx);
  }
  //dihedral_alphabeta1=sqrt(dihedral_alphabeta1/8);
  //dihedral_alphabeta2=sqrt(dihedral_alphabeta2/9);
  */
  //if(noverlaps==0){
  //print(counter, allpositions, r2, rg, dihedral_alphabeta1, dihedral_alphabeta2, allenergies, pofs0, alphahelix_ab, betastrand_ab);
  //}
  //cout << "a" << dihedral_alphabeta1 <<" "<< alphabeta[0][0]  << "\n";
  //cout << "b" << dihedral_alphabeta2 <<" "<< alphabeta[0][1]  << "\n";
  //cout << "hello 7\n";
  //cout << r2 <<" "<< rg<<" r2\n";
}


void initialize(string biasfname, string ab_helix_filename, string ab_sheet_filename, string ab_helix_distance_filename, string ab_sheet_distance_filename,string angle_distance_filename){
  Two_d_grid gridtemplate(gridx, gridy, 2*3.141592654, 2*3.141592654, -3.141592654, -3.141592654);
  //Two_d_grid gridtemplateb(gridx, gridy, 8, 10, 0, 0);
//Nbiases is the number of 2d_grids you want to make, gridx, and gridy are the number of gridpoints you are using
  for(int i=0; i<Nbiases; i++){
    //cout<< i<< " " << Nbiases << "\n";
    bias_grid.push_back(gridtemplate);
    bias_grid[i].readfiles(i+1,biasfname);
    bias_grid[i].setindices();
  }
  //cout << "hello 1\n";
  alphabetabias_helix.readfiles(1,ab_helix_filename);
  alphabetabias_helix.setindices();
  
  alphabetabias_sheet.readfiles(1,ab_sheet_filename);
  alphabetabias_sheet.setindices();
  
  alphabeta_helix_distance.readfiles(1,ab_helix_distance_filename);
  alphabeta_helix_distance.setindices();
  
  alphabeta_sheet_distance.readfiles(1,ab_sheet_distance_filename);
  alphabeta_sheet_distance.setindices();
  
  angle_distance.readfiles(1,angle_distance_filename);
  angle_distance.setindices();
  
  //cout << "index set\n";
  //Nbackbone=3*N_cvs/2 +5;
  Nbackbone=21;
  Nreplicas=1;
  atomistic_configurations.resize(Nreplicas);
  
  for(int i=0; i<Nreplicas; i++){
    atomistic_configurations[i].resize(N_proteins);
    for(int j=0; j<N_proteins; j++){
      atomistic_configurations[i][j].resize(Nbackbone);
      for(int k=0; k<Nbackbone; k++){
	atomistic_configurations[i][j][k].resize(3);
	//cout << atomistic_configurations[i][j][k].size() << "\n";
      }
    }
  }

  beta[0]=0.3400;
  beta[1]=0.3306;
  beta[2]=0.3089;
  beta[3]=0.2887;
  beta[4]=0.2698;
  beta[5]=0.2552;
  beta[6]=0.2357;
  beta[7]=0.2202;
  ofstream colvarfile;
  //colvarfile.open("colvar_pt5rep.data", ios_base::app);
  colvarfile.open(outfile, ios_base::app);
  colvarfile <<"#! FIELDS time"<<"\n";
  colvarfile.close();
  /*crystal_dihedrals[0]=2.41;
  crystal_dihedrals[1]=-2.52;
  crystal_dihedrals[2]=1.98;
  crystal_dihedrals[3]=-1.58;
  crystal_dihedrals[4]=1.98;
  crystal_dihedrals[5]=-0.26;
  crystal_dihedrals[6]=-0.89011;
  crystal_dihedrals[7]=-1.1519;
  crystal_dihedrals[8]=-1.5882;
  crystal_dihedrals[9]=-0.36651;
  crystal_dihedrals[10]=1.4486;
  crystal_dihedrals[11]=0.078;
  crystal_dihedrals[12]=-1.186;
  crystal_dihedrals[13]=1.867;
  crystal_dihedrals[14]=-1.308;
  crystal_dihedrals[15]=2.356;
  crystal_dihedrals[16]=-2.635;*/
  
  crystal_dihedrals[0]=2.41;
  crystal_dihedrals[1]=-2.0;
  crystal_dihedrals[2]=2.3;
  crystal_dihedrals[3]=-1.58;
  crystal_dihedrals[4]=1.98;
  crystal_dihedrals[5]=-0.26;
  crystal_dihedrals[6]=-1.3;
  crystal_dihedrals[7]=-0.8;
  crystal_dihedrals[8]=-1.58;
  crystal_dihedrals[9]=-0.36;
  crystal_dihedrals[10]=1.44;
  crystal_dihedrals[11]=0.5;
  crystal_dihedrals[12]=-1.5;
  crystal_dihedrals[13]=2.2;
  crystal_dihedrals[14]=-1.308;
  crystal_dihedrals[15]=2.356;
  crystal_dihedrals[16]=-2.0;


  alphahelix_dihedrals[0]=-1.047;
  alphahelix_dihedrals[1]=-0.785;
  alphahelix_dihedrals[2]=-1.047;
  alphahelix_dihedrals[3]=-0.785;
  alphahelix_dihedrals[4]=-1.047;
  alphahelix_dihedrals[5]=-0.785;
  alphahelix_dihedrals[6]=-1.047;
  alphahelix_dihedrals[7]=-0.785;
  alphahelix_dihedrals[8]=-1.047;
  alphahelix_dihedrals[9]=-0.785;
  alphahelix_dihedrals[10]=-1.047;
  alphahelix_dihedrals[11]=-0.785;
  alphahelix_dihedrals[12]=-1.047;



  betasheet_dihedrals[0]=-2.44;
  betasheet_dihedrals[1]=2.35;
  betasheet_dihedrals[2]=-2.44;
  betasheet_dihedrals[3]=2.35;
  betasheet_dihedrals[4]=-2.44;
  betasheet_dihedrals[5]=2.35;
  betasheet_dihedrals[6]=-2.44;
  betasheet_dihedrals[7]=2.35;
  betasheet_dihedrals[8]=-2.44;
  betasheet_dihedrals[9]=2.35;
  betasheet_dihedrals[10]=-2.44;
  betasheet_dihedrals[11]=2.35;
  betasheet_dihedrals[12]=-2.44;




}

void print(int counter, vector < vector <double> > & allpositions, vector < vector <double> > & allorientations, vector <double> & alphabeta_ah, vector <double> & alphabeta_bs){
  double distance=0;
  double angle=0;
  ofstream colvarfile;

  colvarfile.open(outfile, ios_base::app);


  colvarfile << counter<<" ";
  for(int i=0; i<N_proteins; i++){
    colvarfile << alphabeta_ah[i];
    colvarfile << " ";
  }

  for(int i=0; i<N_proteins; i++){
    colvarfile << alphabeta_bs[i];
    colvarfile << " ";
  }
  distance=sqrt(pow(allpositions[0][0]-allpositions[1][0],2)+pow(allpositions[0][1]-allpositions[1][1],2)+pow(allpositions[0][2]-allpositions[1][2],2));
  angle = acos(allorientations[0][0]*allorientations[1][0] + allorientations[0][1]*allorientations[1][1] + allorientations[0][2]*allorientations[1][2]);
  colvarfile << distance;
  colvarfile << " ";
  colvarfile << angle;
  colvarfile << " \n";

  colvarfile.close();
  
  //cout << allorientations[0][0]  << " " << allorientations[0][1] << " " << allorientations[0][2] << "\n";

}


