#include <Two_d_grid.h> //includes the 2d grid class for handling everything you want to do with a 2d grid
#include <random>
#include <chrono>
#include <armadillo>
#include <bitset>

//g++ -O2 mc_dihedrals_pt.cpp -I . -I ~/code/armadillo-4.450.2/include -DARMA_USE_BLAS -DARMA_USE_LAPACK -DARMA_DONT_USE_WRAPPER -std=c++0x -lblas -llapack -o mc_dihedrals_pt

using namespace arma;
using namespace std;
int Nbiases=0;
int gridx=200;
int gridy=200;
int N_cvs=0;
double step_size=0.1;
double energy_ref=0.0;
vector <double> beta(8);
vector <double> crystal_dihedrals(50);
vector <vector <double> > alphabeta;
vector <vector <double> > initial_positions;
int Nbackbone=0;
int Nreplicas=0;
int Nsteps=500;
int Nsweeps=1000;
int N_Configs=0;
double kappa=0;
double bessel=0;
vector <vector <int>> binary_configs;
//int Nsteps=0;
//int Nsweeps=0;
vector <Two_d_grid> bias_grid;
Two_d_grid alphabetabias(200,200,8,9,0,0);
void initialize1(string, string);
void initialize2();
void run_mc(int);
void gen_initial(int);
void print(int, vector<double> &, double, double, double, double,  vector<double> &, double);
void backbone(int, vector<double> &, vector<double> &, double);
string outfile;

int main(int argc, char *argv[]){
  string biasfilename;
  string alphabetabiasfname;
  if(argc>1){
    N_cvs = atoi( argv[1]);
    Nbiases = atoi (argv[2]);
    biasfilename = argv[3];
    alphabetabiasfname = argv[4];
    outfile=argv[5];
    kappa = atof(argv[6]);
    bessel = atof(argv[7]);
  }

  
  gen_initial(N_cvs);
  initialize1(biasfilename,alphabetabiasfname);
  for(int i1=0; i1<N_Configs; i1++){//change upper limit back to N_Configs
    initialize2();
    run_mc(i1);
  }

  //run_mc(19114);
  return 0; 
}


void gen_initial(int n_dihedrals){
  int n_configs=0;
  n_configs=pow(2,n_dihedrals);
  binary_configs.resize(n_configs);
  for(int i=0; i<binary_configs.size(); i++){
    binary_configs[i].resize(n_dihedrals);
  }
  
  for(int i=0; i<n_configs; i++){
    bitset<32> A = i;
    for(int j=0; j<n_dihedrals; j++){
      binary_configs[i][j]=A[n_dihedrals-1-j];
    }
  }

  
  binary_configs[0][0]=1;
  binary_configs[0][1]=0;
  binary_configs[0][2]=1;
  binary_configs[0][3]=0;
  binary_configs[0][4]=1;
  binary_configs[0][5]=0;
  binary_configs[0][6]=0;
  binary_configs[0][7]=0;
  binary_configs[0][8]=0;
  binary_configs[0][9]=0;
  binary_configs[0][10]=1;
  binary_configs[0][11]=0;
  binary_configs[0][12]=0;
  binary_configs[0][13]=1;
  binary_configs[0][14]=0;
  binary_configs[0][15]=1;
  binary_configs[0][16]=0;


  N_Configs=n_configs;
}


void run_mc(int config_number){
  vector <vector <double> > positions;
  double pos_temp=0;
  double energy_delta=0;
  double alphabeta1new=0;
  double alphabeta2new=0;
  vector <double> energy(Nreplicas,0);
  vector <double> initialPos(N_cvs);
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
  double pofsold=0;
  double pofsnew=0;
  double pofs0=0;
  double attempts=0;
  double rmax=0;
  double rmin=0;
  double increment=0;
  double rx=0;
  double pofsinc=0;
  double counter=0;
  double partition_function=0;
  double energyav=0;
  double energymin=0;
  double energyps=0;
  double energypsav=0;
  double energyab=0;
  double energyabav=0;
  double ab1av=0;
  double ab2av=0;
  positions.resize(Nreplicas);
  
  alphabeta[0][0]=0;
  alphabeta[0][1]=0;
  for(int k=0; k<positions.size(); k++){
    positions[k].resize(N_cvs);
    for(int j=0; j<N_cvs; j++){
      positions[k][j]=initial_positions[j][binary_configs[config_number][j]];
      initialPos[j]=positions[k][j];
      rx=positions[k][j]-crystal_dihedrals[j];
      if(j<8){
	alphabeta[k][0] += (0.5+0.5*cos(rx));
      }
      else{
	alphabeta[k][1] += (0.5+0.5*cos(rx));
      }
      
    }
  }

  for(int k=0; k<Nreplicas; k++){
    energy[k]=0;
    for(int j=0; j<Nbiases; j++){
      energy[k]+=bias_grid[j].getvalue_linearinterpolation(positions[k][j], positions[k][j+1]);
    }
    energy[k]+=alphabetabias.getvalue_linearinterpolation(alphabeta[k][0],alphabeta[k][1]);
    energyab = alphabetabias.getvalue_linearinterpolation(alphabeta[k][0],alphabeta[k][1]);
    pofsnew=1;
    for(int j=0; j<N_cvs; j++){
      pofsnew *= exp(kappa*cos(positions[k][j]-crystal_dihedrals[j]))/bessel;
    }
    pofsnew+=1;
    energy[k]-=log(pofsnew)/beta[0];
    energyps=-log(pofsnew)/beta[0];
  }
  
  if(config_number==0){
    energy_ref=energy[0];
  }
  for(int i=0; i<=Nsweeps; i++){
    //cout << i << "\n"; 
    for(int k=0; k<Nreplicas; k++){//Nreplicas don't forget to change back
      //cout << k << "\n"; 
      for(int j=0; j<=Nsteps; j++){
	//cout << j << "\n"; 
	dihedral_label = distributionA(generator);
	pos_temp = positions[k][dihedral_label] + (real_distribution(generator)-0.5)*step_size;
	//cout << dihedral_label<< "\n";
	if(pos_temp>xmax-0.01){
	  pos_temp=xmin+(pos_temp-xmax);
	}
	else if(pos_temp<xmin){
	  pos_temp=xmax-(xmin-pos_temp)-0.001;
	}
	
	alphabeta1new=0;
	
	pofsnew=1;
	pofsold=1;

	for(int j1=0; j1<N_cvs; j1++){
	  if(j1==dihedral_label){
	    pofsnew *= exp(kappa*cos(pos_temp-crystal_dihedrals[j1]))/bessel;
	    pofsold *= exp(kappa*cos(positions[k][j1]-crystal_dihedrals[j1]))/bessel;
	  }
	  else{
	    pofsnew *= exp(kappa*cos(positions[k][j1]-crystal_dihedrals[j1]))/bessel;
	    pofsold *= exp(kappa*cos(positions[k][j1]-crystal_dihedrals[j1]))/bessel;
	  }
	}
	pofsnew+=1;
	pofsold+=1;
	pofsinc = (log(pofsold)-log(pofsnew))/beta[0];//beta[0] or beta[k]
	//cout << pofsinc << " " << pofsold << " " << pofsnew << " \n";
	for(int j1=0; j1<8; j1++){
	  if(j1==dihedral_label){
	    rx=pos_temp-crystal_dihedrals[dihedral_label];
	    if(rx>3.14159){
	      rx=-3.14159+(rx-3.14159);
	    }
	    if(rx<-3.14159){
	      rx=3.14159-(-3.14159-rx);
	    }
	  }
	  else{
	    rx=positions[k][j1]-crystal_dihedrals[j1];
	    if(rx>3.14159){
	      rx=-3.14159+(rx-3.14159);
	    }
	    if(rx<-3.14159){
	      rx=3.14159-(-3.14159-rx);
	    }
	  }
	  alphabeta1new+=0.5+0.5*cos(rx);
	}
	//alphabeta1new=sqrt(alphabeta1new/8);
	
	alphabeta2new = 0;
	
	for(int j1=8; j1<N_cvs; j1++){
	  if(j1==dihedral_label){
	    rx=pos_temp-crystal_dihedrals[dihedral_label];
	    if(rx>3.14159){
	      rx=-3.14159+(rx-3.14159);
	    }
	    if(rx<-3.14159){
	      rx=3.14159-(-3.14159-rx);
	    }
	  }
	  else{
	    rx=positions[k][j1]-crystal_dihedrals[j1];
	    if(rx>3.14159){
	      rx=-3.14159+(rx-3.14159);
	    }
	    if(rx<-3.14159){
	      rx=3.14159-(-3.14159-rx);
	    }
	  }
	  alphabeta2new+=0.5+0.5*cos(rx);
	}
	//alphabeta2new=sqrt(alphabeta2new/9);
	//cout << " " << alphabeta1new << " "<<alphabeta2new<<" "<<alphabeta[k][0]<<" "<<alphabeta[k][1]<<" "<< "\n";
	energy_delta=0;

	if(dihedral_label>0 && dihedral_label< N_cvs-1){
	  energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(pos_temp, positions[k][dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(positions[k][dihedral_label],positions[k][dihedral_label+1]);
	  energy_delta += bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[k][dihedral_label-1], pos_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[k][dihedral_label-1],positions[k][dihedral_label]);
	  //cout<< "inc1\n";
	  increment = alphabetabias.getvalue_linearinterpolation(alphabeta1new,alphabeta2new)-alphabetabias.getvalue_linearinterpolation(alphabeta[k][0],alphabeta[k][1]);
	  //cout<<"med " << increment << " " << alphabeta1new << " "<<alphabeta2new<<" "<<alphabeta[k][0]<<" "<<alphabeta[k][1]<<" "<<increment<<"\n";
	  energy_delta += increment + pofsinc;
	}
	else if(dihedral_label==0){	
	  energy_delta = bias_grid[dihedral_label].getvalue_linearinterpolation(pos_temp, positions[k][dihedral_label+1])-bias_grid[dihedral_label].getvalue_linearinterpolation(positions[k][dihedral_label],positions[k][dihedral_label+1]);
	  //cout<< "inc2\n";
	  increment = alphabetabias.getvalue_linearinterpolation(alphabeta1new,alphabeta2new)-alphabetabias.getvalue_linearinterpolation(alphabeta[k][0],alphabeta[k][1]);
	  //cout<<"beg " << increment << " " << alphabeta1new << " "<<alphabeta2new<<" "<<alphabeta[k][0]<<" "<<alphabeta[k][1]<<" "<<increment<<"\n";
	  energy_delta += increment+pofsinc;
	}
	else{
	  energy_delta = bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[k][dihedral_label-1], pos_temp)-bias_grid[dihedral_label-1].getvalue_linearinterpolation(positions[k][dihedral_label-1],positions[k][dihedral_label]);
	  //cout<< "inc2\n";
	  increment = alphabetabias.getvalue_linearinterpolation(alphabeta1new,alphabeta2new)-alphabetabias.getvalue_linearinterpolation(alphabeta[k][0],alphabeta[k][1]);
	  energy_delta += increment+pofsinc;
	  //cout<<"end " << increment << " " << alphabeta1new << " "<<alphabeta2new<<" "<<alphabeta[k][0]<<" "<<alphabeta[k][1]<<" "<<increment<< "\n";
	}
	
	rmax=cos(pos_temp-initialPos[dihedral_label]);
	if((energy_delta<0 && alphabeta1new<7.97 && rmax>0.4) || (real_distribution(generator)<=exp(-beta[k]*energy_delta) && alphabeta1new<7.97 && rmax>0.4)){
	  energy[k] += energy_delta;
	  positions[k][dihedral_label]=pos_temp;
	  alphabeta[k][0]=alphabeta1new;
	  alphabeta[k][1]=alphabeta2new;
	  energyps += pofsinc;
	  energyab += increment;
	  if(k==0){
	    pofs0=pofsnew;
	  }
	}
	else{
	  //cout << "rejected\n";
	  pos_temp=pos_temp;
	  if(k==0){
	    pofs0=pofsold;
	  }
	}
      }
    }
    if(i==0){
      energymin=energy[0];
    }
    if(i%1==0){
      //cout << "hello\n";
      //backbone(i,positions[0],energy,pofs0);
      partition_function+=exp(-(energy[0]-energy_ref)*beta[0]);
      energyav+=energy[0];
      energypsav+=energyps;
      energyabav+=energyab;
      ab1av+=alphabeta[0][0];
      ab2av+=alphabeta[0][1];
      if(energy[0]<energymin){
	energymin=energy[0];
      }
      counter+=1;
      //cout << alphabeta[0][0] << " " << alphabeta[0][1]<<" "<< energy[0] <<" "<<energyps<<" " << energyab << "\n";
    }    
  }
  
  partition_function = -log(partition_function);
  energyav=energyav/counter;
  energypsav=energypsav/counter;
  energyabav=energyabav/counter;
  ab1av=ab1av/counter;
  ab2av=ab2av/counter;
  //cout << config_number << " " << partition_function << " " << energyav << "\n";
  ofstream colvarfile;
  colvarfile.open(outfile, ios_base::app);
  colvarfile << config_number <<" " << partition_function << " " << energyav << " " << energymin << " " << energypsav << " " << energyabav << " "<< ab1av << " " << ab2av << "\n";
  colvarfile.close();
}


void backbone(int counter, vector<double> & allpositions, vector<double> & allenergies, double pofs0){
  vector< vector <double> > backbone;
  vector< double > bond_distances(Nbackbone+8);
  vector< double > bond_angles(Nbackbone+8);
  vector< double > dihedral_angles(Nbackbone+8);
  ofstream backbonefile;
  double dihedral_alphabeta1=0;
  double dihedral_alphabeta2=0;
  string bbstring ("backbone_atoms.xyz");
  backbonefile.open(bbstring,ios::out);
  backbone.resize(Nbackbone);
  for(int k=0; k<Nbackbone; k++){
    backbone[k].resize(3);
  }
  /* backbone[0][0]=0;backbone[0][1]=0;backbone[0][2]=0;
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
  dihedral_angles[25]=2.9367;*/

  /*for(int k=0; k<=25; k++){
    cout << dihedral_angles[k] << "\n";
    }*/

  /*
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
  backbone[2][2]=transformation_matrix(2,3);*/

  //backbonefile <<"N "<<10*backbone[0][0] <<" "<<10*backbone[0][1]<<" "<<10*backbone[0][2]<<"\n";
  //backbonefile <<"C "<<10*backbone[1][0] <<" "<<10*backbone[1][1]<<" "<<10*backbone[1][2]<<"\n";
  //backbonefile <<"C "<<10*backbone[2][0] <<" "<<10*backbone[2][1]<<" "<<10*backbone[2][2]<<"\n";
  /*double cos_ba=0; double sin_ba=0;
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

    //if(k%3==0){
      //backbonefile <<"N " <<10*backbone[k][0] <<" "<<10*backbone[k][1]<<" "<<10*backbone[k][2]<<"\n";
    //}
    //else{
      //backbonefile <<"C " <<10*backbone[k][0] <<" "<<10*backbone[k][1]<<" "<<10*backbone[k][2]<<"\n";
      //}
  
  }*/


  /*double rx, ry, rz=0;
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
  rg = rg/10;*/
  
  double rx=0;
  dihedral_alphabeta2=0; dihedral_alphabeta1=0;
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
  }
  //dihedral_alphabeta1=sqrt(dihedral_alphabeta1/8);
  //dihedral_alphabeta2=sqrt(dihedral_alphabeta2/9);
  double r2=0; double rg=0;
  //print(counter, allpositions, r2, rg, dihedral_alphabeta1, dihedral_alphabeta2, allenergies, pofs0);
  //cout << "a" << dihedral_alphabeta1 <<" "<< alphabeta[0][0]  << "\n";
  //cout << "b" << dihedral_alphabeta2 <<" "<< alphabeta[0][1]  << "\n";
  //cout << "hello 7\n";
  //cout << r2 <<" "<< rg<<" r2\n";
}


void initialize1(string biasfname, string alphabetabiasfname){
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
  //cout << "hello 1.5\n";
  alphabetabias.setindices();
  //cout << "index set\n";
  Nbackbone=3*N_cvs/2 +5;
  Nreplicas=1;
  beta[0]=0.353742;
  //ofstream colvarfile;
  //colvarfile.open("colvar_pt5rep.data", ios_base::app);
  //colvarfile.open(outfile, ios_base::app);
  //colvarfile <<"#! FIELDS time psi-1 phi-2 psi-2 phi-3 psi-3 psi-4 phi-5 psi-5 phi-6 psi-6 phi-7 psi-7 phi-8 psi-8 phi-9 psi-9 phi-10 d1 rg dihedral_alphabeta1 dihedral_alphabeta2 bias"<<"\n";
  //colvarfile.close();
  crystal_dihedrals[0]=2.41;
  crystal_dihedrals[1]=-2.0;
  crystal_dihedrals[2]=2.3;
  crystal_dihedrals[3]=-1.58;
  crystal_dihedrals[4]=1.98;
  // crystal_dihedrals[5]=-1.34;
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

  initial_positions.resize(N_cvs);
  for(int k=0; k<N_cvs; k++){
    initial_positions[k].resize(2);
  }

  initial_positions[0][0]=-1.5;
  initial_positions[0][1]=3.1;
  initial_positions[1][0]=-1.5;
  initial_positions[1][1]=1.2;
  initial_positions[2][0]=-0.3;
  initial_positions[2][1]=2.5;
  initial_positions[3][0]=-1.5;
  initial_positions[3][1]=1.15;
  initial_positions[4][0]=-0.7;
  initial_positions[4][1]=2.4;
  initial_positions[5][0]=-0.3;
  initial_positions[5][1]=2.4;
  initial_positions[6][0]=-1.4;
  initial_positions[6][1]=1.1;
  initial_positions[7][0]=-0.2;
  initial_positions[7][1]=2.4;
  initial_positions[8][0]=-1.5;
  initial_positions[8][1]=1.1;
  initial_positions[9][0]=-0.3;
  initial_positions[9][1]=2.1;
  initial_positions[10][0]=-1.5;
  initial_positions[10][1]=1.3;
  initial_positions[11][0]=0.0;
  initial_positions[11][1]=3.0;
  initial_positions[12][0]=-1.5;
  initial_positions[12][1]=1.2;
  initial_positions[13][0]=-0.3;
  initial_positions[13][1]=2.2;
  initial_positions[14][0]=-1.4;
  initial_positions[14][1]=1.2;
  initial_positions[15][0]=-0.3;
  initial_positions[15][1]=2.3;
  initial_positions[16][0]=-1.4;
  initial_positions[16][1]=1.1;

  alphabeta.resize(Nreplicas);
  for(int k=0; k<alphabeta.size(); k++){
    alphabeta[k].resize(2);
  }

}

void initialize2(){
  double rx=0;

}

void print(int counter, vector<double> & allpositions, double r2, double rg, double dihedral_rmsd1, double dihedral_rmsd2, vector<double> & allenergies, double pofs0){
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
  colvarfile << pofs0 << " ";
  colvarfile << "\n";

  colvarfile.close();

}


