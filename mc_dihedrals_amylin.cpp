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
double pos_step_size=0.2;
double orientation_step_size=0.3;
vector <double> beta(8);
vector <double> crystal_dihedrals(50);
vector <double> alphahelix_dihedrals(50);
vector <double> betastrand_dihedrals(50);
vector <vector <double> > alphabeta;
vector <double> initial_dihedrals(50);
int Nbackbone=0;
int Nreplicas=0;
int Nsteps=8000;
int Nsweeps=5000000;
int N_proteins=0;
double kappa=0;
double bessel=0;
double fold_weight=0;
//int Nsteps=0;
//int Nsweeps=0;

vector <Two_d_grid> bias_grid;
Two_d_grid alphabetabias(200,200,8,9,0,0);
void initialize(string, string);
void run_mc();
void print(int, vector<double> &, double, double, double, double,  vector<double> &, double, double, double);
void backbone(int, vector<double> &, vector<double> &, double);
string outfile;

int main(int argc, char *argv[]){
  string biasfilename;
  string alphabetabiasfname;
  if(argc>1){
    N_cvs = atoi( argv[1]);
    Nbiases = atoi (argv[2]);
    biasfilename = argv[3];
    N_proteins = atoi (argv[4]);
    alphabetabiasfname = argv[4];
    outfile=argv[5];
    kappa = atof(argv[6]);
    bessel = atof(argv[7]);
    fold_weight = atof(argv[8]);
  }
  cout << "hello\n";
  initialize(biasfilename,alphabetabiasfname);
  cout << "hello\n";
  run_mc();
  
  return 0; 
}

void run_mc(){
  vector <vector <vector <double> > > dihedral_angles(Nreplicas, vector <vector <double> >(N_proteins, vector <double>(N_cvs))); //index order: 1st corresponds to replica, 2nd to protein, 3rd to dihedral 
  vector <vector <vector <double> > > positions(Nreplicas, vector <vector <double> >(N_proteins, vector <double>(3)));
  vector <vector <vector <double> > > orientations(Nreplicas, vector <vector <double> >(N_proteins, vector <double>(3)));
  vector <vector <double> > ab_alphahelix(Nreplicas, vector <double >(N_proteins));
  vector <vector <double> > ab_betasheet(Nreplicas, vector <double >(N_proteins));

  vector<double> pos_temp(3);
  vector<double> orientation_temp(3);
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
  double pofsnew=0;
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
      
      positions[k][j][0]=2.07*real_distribution(generator);
      positions[k][j][1]=2.07*real_distribution(generator);
      positions[k][j][2]=2.07*real_distribution(generator);
      
      orientations[k][j][0]=0;
      orientations[k][j][1]=0;
      if(j%2==0){
	orientations[k][j][2]=1;
      }
      else{
	orientations[k][j][2]=-1;
      }
      ab_alphahelix[k][j]=0;
      ab_betasheet[k][j]=0;
      for(int j1=0; j1<dihedral_angles[k][j].size(); j1++){
	dihedral_angles[k][j][j1]=-3.14+real_distribution(generator)*2*3.14;
	double angle_delta = dihedral_angles[k][j][j1]-alphahelix_dihedrals[j1];
	ab_alphahelix[k][j] += 0.5+0.5*cos(angle_delta);
	angle_delta = dihedral_angles[k][j][j1]-betasheet_dihedrals[j1];
	ab_betasheet[k][j] += 0.5+0.5*cos(angle_delta);
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
	
	alphabeta1new=0;
	
	pofsnew=1;
	pofsold=1;

	for(int j1=0; j1<N_cvs; j1++){
	  if(j1==dihedral_label){
	    pofsnew *= exp(kappa*cos(angle_temp-crystal_dihedrals[j1]))/bessel;
	    pofsold *= exp(kappa*cos(dihedral_angles[k][j1]-crystal_dihedrals[j1]))/bessel;
	  }
	  else{
	    pofsnew *= exp(kappa*cos(dihedral_angles[k][j1]-crystal_dihedrals[j1]))/bessel;
	    pofsold *= exp(kappa*cos(dihedral_angles[k][j1]-crystal_dihedrals[j1]))/bessel;
	  }
	}
	
	pofsnew *= fold_weight;
	pofsold *= fold_weight;
	
	pofsnew+=1;
	pofsold+=1;
	pofsinc = (log(pofsold)-log(pofsnew))/beta[0];//beta[0] or beta[k]
	//cout << pofsinc << " " << pofsold << " " << pofsnew << " \n";
	for(int j1=0; j1<12; j1++){
	  if(j1==dihedral_label){
	    rx=angle_temp-alphahelix_dihedrals[dihedral_label];
	    rx2=angle_temp-betasheet_dihedrals[dihedral_label];
	    if(rx>3.14159){
	      rx=-3.14159+(rx-3.14159);
	    }
	    if(rx<-3.14159){
	      rx=3.14159-(-3.14159-rx);
	    }
	  }
	  else{
	    rx=dihedral_angles[k][protein_label][j1]-alphahelix_dihedrals[j1];
	    rx2=dihedral_angles[k][protein_label][j1]-betasheet_dihedrals[j1];
	    if(rx>3.14159){
	      rx=-3.14159+(rx-3.14159);
	    }
	    if(rx<-3.14159){
	      rx=3.14159-(-3.14159-rx);
	    }
	  }
	  ab_alphahelix_new+=0.5+0.5*cos(rx);
	  ab_betasheet_new+=0.5+0.5*cos(rx2);
	}

	energy_delta=0;
	
	if(dihedral_label>0 && dihedral_label< N_cvs-1){
	  energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(angle_temp, dihedral_angles[k][protein_label][dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label],dihedral_angles[k][protein_label][dihedral_label+1]);
	  energy_delta += bias_grid[dihedral_label-1].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label-1], angle_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label-1],dihedral_angles[k][dihedral_label]);

	}
	else if(dihedral_label==0){	
	  energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(angle_temp, dihedral_angles[k][protein_label][dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label],dihedral_angles[k][protein_label][dihedral_label+1]);

	}
	else{
	  energy_delta = bias_grid[dihedral_label-1].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label-1], angle_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(dihedral_angles[k][protein_label][dihedral_label-1],dihedral_angles[k][protein_label][dihedral_label]);

	}
	
	energy_delta += alphabetabias_helix.getvalue_linearinterpolation(ab_alphahelix_new)-alphabetabias_helix.getvalue_linearinterpolation(ab_alphahelix[k][protein_label]);

	energy_delta += alphabetabias_sheet.getvalue_linearinterpolation(ab_betasheet_new)-alphabetabias_sheet.getvalue_linearinterpolation(ab_betasheet[k][protein_label]);

	for(int i1=0; i1< N_proteins; i1++){
	  if(i1 != protein_label){
	    double ip_r=0;
	    double ip_r_new=0;
	    double interproteinx = positions[k][protein_label][0]-positions[k][i1][0];
	    double interproteiny = positions[k][protein_label][1]-positions[k][i1][1];
	    double interproteinz = positions[k][protein_label][2]-positions[k][i1][2];
	    double interproteinx_new = postemp[0]-positions[k][i1][0];
	    double interproteiny_new = postemp[1]-positions[k][i1][1];
	    double interproteinz_new = postemp[2]-positions[k][i1][2];
	    ip_r = pow(interproteinx*interproteinx + interproteiny*interproteiny + interproteinz*interproteinz, 0.5);
	    ip_r_new = pow(interproteinx_new*interproteinx_new + interproteiny_new*interproteiny_new + interproteinz_new*interproteinz_new, 0.5);


	  }
	  

	  
	  energy_delta += ;

	}

	energy_delta += pofsinc;


	if((energy_delta<0 && alphabeta1new<7.98 && alphabeta2new<8.98) || (real_distribution(generator)<=exp(-beta[k]*energy_delta) && alphabeta1new<7.98 && alphabeta2new<8.98)){
	  energy[k] += energy_delta;
	  dihedral_angles[k][protein_label][dihedral_label]=angle_temp;
	  positions[k][protein_label][0]=pos_temp[0];
	  positions[k][protein_label][1]=pos_temp[1];
	  positions[k][protein_label][2]=pos_temp[2];
	  
	  orientations[k][protein_label][0]=orientation_temp[0];
	  orientations[k][protein_label][1]=orientation_temp[1];
	  orientations[k][protein_label][2]=orientation_temp[2];
	  
	  
	  ab_alphahelix[k][protein_label]=ab_alphahelix_new;
	  ab_betasheet[k][protein_label]=ab_betasheet_new;

	  alphabeta[k][0]=alphabeta1new;

	  alphabeta[k][1]=alphabeta2new;
	  if(k==0){
	    pofs0=pofsnew;
	  }
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
      //cout << "hello\n";
      backbone(i,dihedral_angles[0],energy,pofs0);
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

void backbone(int counter, vector<double> & allpositions, vector<double> & allenergies, double pofs0){
  vector< vector <double> > backbone;
  vector< double > bond_distances(Nbackbone+8);
  vector< double > bond_angles(Nbackbone+8);
  vector< double > dihedral_angles(Nbackbone+8);
  ofstream backbonefile;
  double dihedral_alphabeta1=0;
  double dihedral_alphabeta2=0;
  double alphahelix_ab=0;
  double betastrand_ab=0;
  string bbstring ("backbone_atoms.xyz");
  backbonefile.open(bbstring,ios::out);
  backbone.resize(Nbackbone);
  for(int k=0; k<Nbackbone; k++){
    backbone[k].resize(3);
  }
  backbone[0][0]=0;backbone[0][1]=0;backbone[0][2]=0;
  bond_distances[0]=0.133386;
  for(int k=0; k<=N_cvs/2+2; k++){
    bond_distances[3*k]=0.147101;
    bond_distances[1+3*k]=0.1546;
    bond_distances[2+3*k]=0.133649;
  }

  bond_distances[0]=0.1606;
  bond_distances[1]=0.1525;
  bond_distances[2]=0.1399;
  bond_distances[3]=0.14909;
  bond_distances[4]=0.14855;
  bond_distances[5]=0.129112;
  bond_distances[6]=0.144655;
  bond_distances[7]=0.148;
  bond_distances[8]=0.1327;
  bond_distances[9]=0.1456;
  bond_distances[10]=0.1515;
  bond_distances[11]=0.1328;
  bond_distances[12]=0.1460;
  bond_distances[13]=0.1529;
  bond_distances[14]=0.1333;
  bond_distances[15]=0.1426;
  bond_distances[16]=0.1483;
  bond_distances[17]=0.1419;
  bond_distances[18]=0.1399;
  bond_distances[19]=0.1417;
  bond_distances[20]=0.1288;
  bond_distances[21]=0.1495;
  bond_distances[22]=0.1503;
  bond_distances[23]=0.1333;
  bond_distances[24]=0.1442;
  bond_distances[25]=0.1536;
  bond_distances[26]=0.13601;
  bond_distances[27]=0.1458;
  bond_distances[28]=0.15229;

  bond_angles[0]=2.1957;
  for(int k=0; k<=N_cvs/2+2; k++){
    bond_angles[3*k]=2.0031;
    bond_angles[1+3*k]=2.0653;
    bond_angles[2+3*k]=2.17011;
  }
  bond_angles[0]=1.8264;
  bond_angles[1]=2.012;
  bond_angles[2]=2.171;
  bond_angles[3]=1.8646;
  bond_angles[4]=2.1627;
  bond_angles[5]=2.114;
  bond_angles[6]=1.9621;
  bond_angles[7]=2.1412;
  bond_angles[8]=2.214;
  bond_angles[9]=1.8408;
  bond_angles[10]=2.163;
  bond_angles[11]=2.328;
  bond_angles[12]=1.9584;
  bond_angles[13]=2.018;
  bond_angles[14]=2.133;
  bond_angles[15]=2.0805;
  bond_angles[16]=1.9535;
  bond_angles[17]=2.115;
  bond_angles[18]=2.088;
  bond_angles[19]=2.087;
  bond_angles[20]=2.2922;
  bond_angles[21]=1.908;
  bond_angles[22]=2.0831;
  bond_angles[23]=2.1431;
  bond_angles[24]=1.875;
  bond_angles[25]=2.1307;
  bond_angles[26]=2.1735;
  bond_angles[27]=2.007;



  int j1=0;

  for(int k=0; k<Nbackbone/3; k++){
    dihedral_angles[3*k+0]=allpositions[j1];j1++;
    dihedral_angles[3*k+1]=3.14159265;
    if(k==2){
      dihedral_angles[3*k+2]=-1.34;
    }
    else{
      dihedral_angles[3*k+2]=allpositions[j1];j1++;
    }
  }

  dihedral_angles[1]=2.971;
  dihedral_angles[4]=2.961;
  dihedral_angles[7]=-2.849;
  dihedral_angles[10]=-2.89;
  dihedral_angles[13]=-3.1265;
  dihedral_angles[16]=-3.0913;
  dihedral_angles[19]=-3.1074;
  dihedral_angles[22]=2.9107;
  dihedral_angles[25]=2.9367;

  /*for(int k=0; k<=25; k++){
    cout << dihedral_angles[k] << "\n";
    }*/

  
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

    backbone[k][0]=transformation_matrix(0,3);
    backbone[k][1]=transformation_matrix(1,3);
    backbone[k][2]=transformation_matrix(2,3);

    /*if(k%3==0){
      backbonefile <<"N " <<10*backbone[k][0] <<" "<<10*backbone[k][1]<<" "<<10*backbone[k][2]<<"\n";
    }
    else{
      backbonefile <<"C " <<10*backbone[k][0] <<" "<<10*backbone[k][1]<<" "<<10*backbone[k][2]<<"\n";
      }*/
  
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
  if(noverlaps==0){
    print(counter, allpositions, r2, rg, dihedral_alphabeta1, dihedral_alphabeta2, allenergies, pofs0, alphahelix_ab, betastrand_ab);
  }
  //cout << "a" << dihedral_alphabeta1 <<" "<< alphabeta[0][0]  << "\n";
  //cout << "b" << dihedral_alphabeta2 <<" "<< alphabeta[0][1]  << "\n";
  //cout << "hello 7\n";
  //cout << r2 <<" "<< rg<<" r2\n";
}


void initialize(string biasfname, string alphabetabiasfname){
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
  alphabetabias.readfiles(1,alphabetabiasfname);
  cout << "hello 1.5\n";
  alphabetabias.setindices();
  cout << alphabetabias.getvalue_linearinterpolation(4,5) << "\n";
  //cout << "index set\n";
  Nbackbone=3*N_cvs/2 +5;
  Nreplicas=1;
  beta[0]=0.353742;
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
  colvarfile <<"#! FIELDS time psi-1 phi-2 psi-2 phi-3 psi-3 psi-4 phi-5 psi-5 phi-6 psi-6 phi-7 psi-7 phi-8 psi-8 phi-9 psi-9 phi-10 d1 rg dihedral_alphabeta1 dihedral_alphabeta2 dihedral_alphahelix dihedral_betastrand bias energy"<<"\n";
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



  betastrand_dihedrals[0]=-2.44;
  betastrand_dihedrals[1]=2.35;
  betastrand_dihedrals[2]=-2.44;
  betastrand_dihedrals[3]=2.35;
  betastrand_dihedrals[4]=-2.44;
  betastrand_dihedrals[5]=2.35;
  betastrand_dihedrals[6]=-2.44;
  betastrand_dihedrals[7]=2.35;
  betastrand_dihedrals[8]=-2.44;
  betastrand_dihedrals[9]=2.35;
  betastrand_dihedrals[10]=-2.44;
  betastrand_dihedrals[11]=2.35;
  betastrand_dihedrals[12]=-2.44;


  initial_dihedrals[0]=-1.5;
  initial_dihedrals[1]=-1.5;
  initial_dihedrals[2]=2.5;
  initial_dihedrals[3]=-1.5;
  initial_dihedrals[4]=-0.7;
  initial_dihedrals[5]=2.4;
  initial_dihedrals[6]=-1.4;
  initial_dihedrals[7]=2.4;
  initial_dihedrals[8]=-1.5;
  initial_dihedrals[9]=2.1;
  initial_dihedrals[10]=-1.5;
  initial_dihedrals[11]=3;
  initial_dihedrals[12]=-1.5;
  initial_dihedrals[13]=2.2;


  alphabeta.resize(Nreplicas);
  for(int k=0; k<alphabeta.size(); k++){
    alphabeta[k].resize(2);
  }
  double rx=0;
  for(int j=0; j<Nreplicas; j++){
    for(int k=0; k<N_cvs; k++){
      rx=0-crystal_dihedrals[k];
      rx=0;//start from crystal structure
      rx=initial_dihedrals[k]-crystal_dihedrals[k];
      /*if(k<15){
	rx=0.3;
	}*/
      if(k<8){
	//dihedral_rmsd1+=rx*rx;
	alphabeta[j][0]+=(0.5+0.5*cos(rx));
	//alphabeta[j][0]+=rx*rx;
      }
      else{
	//dihedral_rmsd2+=rx*rx;
	alphabeta[j][1]+=(0.5+0.5*cos(rx));
	//alphabeta[j][1]+=rx*rx;
      }
    }
    //alphabeta[j][0]=sqrt(alphabeta[j][0]/8);
    //alphabeta[j][1]=sqrt(alphabeta[j][1]/9);
    //alphabeta[j][0]+=1;
    cout << alphabeta[j][0] << " " << alphabeta[j][1]<< "\n";
  }
  

}

void print(int counter, vector<double> & allpositions, double r2, double rg, double dihedral_rmsd1, double dihedral_rmsd2, vector<double> & allenergies, double pofs0, double alphahelix_ab, double betastrand_ab){
  ofstream colvarfile;
  //colvarfile.open("colvar_pt5rep.data", ios_base::app);
  colvarfile.open(outfile, ios_base::app);
  colvarfile << counter<<" ";
  for(int i=0; i<N_cvs; i++){
    colvarfile << allpositions[i];
    colvarfile << " ";
  }
  colvarfile << r2 <<" ";
  colvarfile << rg <<" ";
  colvarfile << dihedral_rmsd1 <<" ";
  colvarfile << dihedral_rmsd2 <<" ";
  colvarfile << alphahelix_ab <<" ";
  colvarfile << betastrand_ab <<" ";
  colvarfile << pofs0 << " ";
  colvarfile << allenergies[0] << " ";
  colvarfile << "\n";

  colvarfile.close();

}


