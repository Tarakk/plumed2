/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/SwitchingFunction.h"
#include "tools/KernelFunctions.h"
#include "tools/Communicator.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>

using namespace std;

namespace PLMD{
  namespace colvar{
    
//+PLUMEDOC COLVAR FEDERQ_CLUST
/*

*/
//+ENDPLUMEDOC
   
    class GsmacShell : public Colvar{
      SwitchingFunction distanceswitch;  // definition of switching function for f_ij
      SwitchingFunction coordswitch;     // definition of switching function for rho_i

    public:

      vector<AtomNumber> start;      // LIST OF MOLECULES NEEDED BY SMACK -> START
      vector<AtomNumber> end;        // LIST OF MOLECULES NEEDED BY SMACK -> END
      vector<AtomNumber> center1;    // LIST OF MOLECULES NEEDED BY SMACK -> CENTER1
      //vector<AtomNumber> center2;    // LIST OF MOLECULES NEEDED BY SMACK -> CENTER2
      vector<AtomNumber> all_atoms;  // LIST OF MOLECULES NEEDED BY SMACK -> CENTER
      vector<double> angles;
      vector<double> width;

      unsigned int mols;   // total number of molecules involved in the calculations
      unsigned int atoms;  // total number of atoms involved in the calculation

      double drbcut, routcr, gammab, sigmaoutcr, f_lc, f_ll; // values for kernel function
      vector<double> kval {0.0, 0.0, 0.0};  // initialize kernel function and its derivative
      vector<double> dkval {0.0, 0.0, 0.0};  // initialize kernel function and its derivative
      void kernel(double);  // kernel function

      ofstream fdbg;   // definition of object needed for debugging
      ofstream rdbg;
      
      
      struct GsmacShelllist_s {   // STRUCTURE NEEDED TO CONSTRUCT THE VERLET LIST
	double rcut;      // CUT OFF OF THE LIST
	double rskin;     // SKIN TO CHECK IF THE LIST HAS TO BE RECALCULATED
	double rskin2;    // SKIN SQUARED
	int step;
	vector<int> nn;       // NUMBERS OF NEIGHBORS
	vector<vector<int> > ni; // List, it's a vector of vectors, for each molecule one vector
	vector<Vector> pos1;  // POSITIONS CENTER MOLECULE
	vector<Vector> pos2;  // POSITIONS NEIGHBORING MOLECULES
      } GsmacShelllist;  // object of type GsmacShelllist_s defined
      
      
      GsmacShell(const ActionOptions&);              //    CONSTRUCTOR
      ~GsmacShell();                                 //    DESTRUCTOR
      // active methods:
      static void registerKeywords( Keywords& keys );  // KEYWORDS
      virtual void calculate();                        // CALCULATE CV 

      void GsmacShell_newlist(vector<AtomNumber>&, GsmacShelllist_s &GsmacShelllist);   // NEW VERLIST
      void GsmacShell_checklist(vector<AtomNumber>&, GsmacShelllist_s &GsmacShelllist); // CHECK THE VERLET LIST
      double dotprod(Vector,Vector);
      double norm2(Vector);
      
    };
    
    PLUMED_REGISTER_ACTION(GsmacShell,"GSMACSHELL")
    
    void GsmacShell::registerKeywords( Keywords& keys ){

      Colvar::registerKeywords(keys);
      keys.add("atoms","CENTER1","the labels of the atoms acting as center of the molecules");
      //keys.add("atoms","CENTER2","the labels of the atoms acting as center of the neighboring molecules");
      keys.add("atoms","START","the labels of the atoms acting as start of the intramolecular vector");
      keys.add("atoms","END","the labels of the atoms acting as end  of the intramolecular vector");
      keys.add("compulsory","ANGLES"," Angles that have to be used in the SMAC calculation (need to be associated with width)");
      keys.add("compulsory","N_ANGLES"," Number of angles that have to be used in the SMAC calculation (need to be associated with width)");
      keys.add("compulsory","WIDTH"," Width of the Gaussian used to calculate the angles distribution (need to be associated with angles)");
      keys.add("optional","SWITCH","Used if you want to employ an alternative swiching function.");
      keys.add("optional","SWITCH_COORD","Used if you want to employ an alternative swiching function.");
      keys.add("compulsory","R_CUT","Verlet lists, default r_cut = R_max");
      keys.add("compulsory","R_SKIN","Verlet lists, default r_skin = 1.25 R_max");
     
      //keys.add("compulsory","ROUTCR","outer boundary of the control region shell");
      //keys.add("compulsory","SIGMAOUTCR","outer CR side width of gaussian for the interval switching function");
      keys.add("compulsory","DRBCUT","up to this region bias will be imposed");
      keys.add("compulsory","GAMMAB","parameter for Fermi switching function used to select particles within drbcut");

      keys.remove("NOPBC");

    }
    
    GsmacShell::GsmacShell(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {
      
      string sw,errors;        //// SWITCH TO DEFINE WHICH SWF WE ARE GOING TO USE FOR THE PAIR
      parse("SWITCH",sw);
      if(sw.length()>0){
	distanceswitch.set(sw,errors);
	if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
      } else {                  // IN CASE NO SF ARE PROVIDED, USE THE STANDARD RATIONAL ONE
	int nn=6;
	int mm=12;
	double d0=0.0;
	double r0=0.0;
	parse("R_0",r0);
	if(r0<=0.0) error("R_0 should be explicitly specified and positive");
	parse("D_0",d0);
	parse("NN",nn);
	parse("MM",mm);
	distanceswitch.set(nn,mm,r0,d0);
      }
      
      parse("SWITCH_COORD",sw);   //// SWITCH TO DEFINE WHICH SWF WE ARE GOING TO USE FOR THE DENSITY
      if(sw.length()>0){
	coordswitch.set(sw,errors);
	if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
      } else {                  // IN CASE NO SF ARE PROVIDED, USE THE STANDARD RATIONAL ONE
	int nn=6;
	int mm=12;
	double d0=0.0;
	double r0=0.0;
	parse("R_0",r0);
	if(r0<=0.0) error("R_0 should be explicitly specified and positive");
	parse("D_0",d0);
	parse("NN",nn);
	parse("MM",mm);
	coordswitch.set(nn,mm,r0,d0);
      }


      unsigned int nAngles;
      parse("N_ANGLES",nAngles);
      angles.resize(nAngles);
      parseVector("ANGLES",angles);
      log.printf("ANGLES ARE: %f %f %f %f \n", angles[0], angles[1], angles[2], angles[3]);
      width.resize(nAngles);
      parseVector("WIDTH",width);

      parseAtomList("CENTER1",center1);
      //parseAtomList("CENTER2",center2);
      parseAtomList("START",start);
      parseAtomList("END",end);

      parse("DRBCUT",drbcut);
      //parse("ROUTCR",routcr);
      parse("GAMMAB",gammab);
      //parse("SIGMAOUTCR",sigmaoutcr);

      mols=center1.size();
      GsmacShelllist.pos1.resize(mols);
      GsmacShelllist.pos2.resize(mols);
      GsmacShelllist.nn.resize(mols);
      GsmacShelllist.ni.resize(mols);
      for(unsigned int i=0; i < mols ; i++){
	GsmacShelllist.ni[i].resize(mols);
      }
      
      
      if(mols==0) error("no molecules specified");
      
      for(unsigned int i=0 ; i < mols ; i++){
	all_atoms.push_back(center1[i]);  // vector.push_back(element) adds element to the
      }                                   // end of the vector
      for(unsigned int i=0 ; i < mols ; i++){
	all_atoms.push_back(center1[i]);    // add the center twice to all_atoms vector
      }
      for(unsigned int i=0 ; i < mols ; i++){
	all_atoms.push_back(start[i]);
      }
      for(unsigned int i=0 ; i < mols ; i++){
	all_atoms.push_back(end[i]);  // all_atoms is now a vector of size 4*mols
      }
      
      atoms=all_atoms.size();
      
      addValueWithDerivatives();  // informs plumed core that we require space to store the value      
      setNotPeriodic();           // of the CV and that the CV will act on the list of atoms
      requestAtoms(all_atoms);    // named all_atoms

      double r_cut,r_skin;
      parse("R_CUT",r_cut);
      parse("R_SKIN",r_skin);

      GsmacShelllist.rcut=r_cut;
      GsmacShelllist.rskin=r_skin;
      GsmacShelllist.rskin2=r_skin*r_skin;
      GsmacShelllist.step=0;
     
      checkRead();  // check that everything on the input line has been read properly,
                    // all parse command should follow before checkRead()
      
      log<<"  contacts are counted with cutoff "<< distanceswitch.description()<<"\n";
      log<<"  correct density is evaluated above with cutoff "<< coordswitch.description()<<"\n";

    }
    
void GsmacShell::kernel(double drb) {
  // switching function to limit action of Gsmac CV


// transform input values into cartesian coordinates
//double zbox = getBox()[2][2];
//double drbcut_c = drbcut*zbox;
//double routcr_c = routcr*zbox;
//double gammab_c = gammab*zbox;
//double sigmaoutcr_c = sigmaoutcr*zbox;

// Instead of z, here radial distance 'd' should be used
  double drbcut_c = drbcut;
  double gammab_c = gammab;

  double f_lc = 0;
  double f_ll = 0;

  // get value of switching function in dependence of z
  // f_ll = 1/(1 + std::exp(-sigmaoutcr_c*(z-routcr_c)));
  f_lc = 1/(1 + std::exp(gammab_c*(drb-drbcut_c)));

  kval[2] = f_lc*(1-f_ll);
  kval[1] = kval[2];
  kval[0] = kval[2];
  dkval[2] = f_lc*(1-f_lc)*(1-f_ll) - f_lc*f_ll*(1-f_ll);

}


    
void GsmacShell::calculate()
{
  
  double cv_val;   // CV
  cv_val=0;
  
  Tensor virial;   // VIRIAL
  vector<Vector> deriv(getNumberOfAtoms());  // DERIVATIVES, vector of customized Plumed vectors
  
  if(GsmacShelllist.step==0) GsmacShell_newlist(center1,GsmacShelllist); // CHECK IF NEIGHBOR LIST HAVE TO BE CONSTRUCTED
  GsmacShelllist.step=1; GsmacShell_checklist(center1,GsmacShelllist);   // CALL NEIGHBOR LIST IN CASE
  
  unsigned int stride=comm.Get_size();  // SET THE PARALLELIZATION VARIABLES for the for loops
  unsigned int rank=comm.Get_rank();
  
  stride=comm.Get_size();
  rank=comm.Get_rank();
  
  // get the vector of each molecule
  vector<Vector> v(getNumberOfAtoms()); 
  for(unsigned int i = rank; i < mols; i += stride) {  // DEFINE THE INTRAMOLECULAR VECTOR
    int start=i + (int) mols + (int) mols;  // (int) is added to mols because mols was declared before as unsigned int
    int end=i + (int) mols + (int) mols + (int) mols;
    v[i]=pbcDistance(getPosition(start),getPosition(end));  // THIS IS DONE ONCE PER STEP FOR ALL MOLECULES
  }
  
  comm.Sum(v);  // MERGING UP

  for(unsigned int i=rank;i<mols;i+=stride) {                   // SUM OVER MOLECULES
    double n,angtot;
    int start_i,end_i;
    double S1,S2,S3;
    start_i=start[i].serial();
    end_i=end[i].serial();

    
    vector<double> f(GsmacShelllist.nn[i]);                            // SWITCHING FUNCTION
    vector<double> omega(GsmacShelllist.nn[i]);                            // SWITCHING FUNCTION
    vector<Vector> domega1(GsmacShelllist.nn[i]);                      // ANGULAR PART1
    vector<Vector> domega2(GsmacShelllist.nn[i]);                      // ANGULAR PART2
    vector<Vector> df(GsmacShelllist.nn[i]);                      // ANGULAR PART2
   

 // identify the atoms within r_cut and calculate their coordination numbers
    n=0.;
    angtot=0.;
    
    Vector dist;   // distance between the pair of atoms 

    for( int j=0; j < GsmacShelllist.nn[i]; ++j) {                   // SUM OVER NEIGHBORS
      int index_j;
      index_j=GsmacShelllist.ni[i][j]-mols;                                      // TRUE INDEX OF THE J NEIGHBOR
      //      if(i==index_j){continue;}
      Vector comp;
      double modij;
      dist=pbcDistance(getPosition(i),getPosition(index_j+mols));    // DISTANCE BETWEEN THEM
      modij=dist.modulo();                                              //
      
      double dfunc=0.;                                                        // CALCULATING SWITCHING FUNCTION AND 
      f[j] = distanceswitch.calculate(modij,dfunc);
      n += f[j];                // CALCULATING THE COORDINATION NUMBER

      
      for(unsigned int ix=0 ;ix<3 ;ix++){ 	
	df[j][ix] = -dfunc * dist[ix];   // DERIVATIVE OF THE SWITCHING FUNCTION
      }
      
      
  // angular contribution to the CV
      // ANGLE BETWEEN THE DIRECTIVE VECTOR        THE INDEX ARE WRONG IN THE SCALAR PRODUCT! CHANGE THEM!                                                   
      double alpha;
      alpha = dotprod(v[i],v[index_j])/sqrt(norm2(v[i])*norm2(v[index_j]));    ///  CHANGE THIS PART IN MOR C++ STYLE
      
      if(alpha>=1.0){                                                                // NOT SURE IF THIS CHECK IS NECESSARY ANYMORE. DOUBLECHECK!
	alpha=0.99999;
      }else if(alpha<=-1.0){
	alpha=-0.99999;
      }

      double phi;
      double domega;
      phi=acos(alpha);
      
      omega[j]=0.;
      domega=0.;
      for(unsigned int k=0 ; k < angles.size() ;k++){        // ANGULAR SERIES! NEED TO BE SUBSTITUTED WITH THE KERNEL
	double e1;	
	e1 = exp(-((phi - angles[k])*(phi - angles[k]))/(2*width[k]*width[k]));    // GAUSSIAN
	omega[j]+=e1;                                                                // SUM OVE GAU
	domega += - e1*(phi - angles[k])/(width[k]*width[k]);                        // DERIV OF GAU
      }
      
      
      double dvac;                                                                  // PART OF THE DERIVATIVES OF THE ANGULAR PART
      dvac = - 1/(sqrt((norm2(v[i])*norm2(v[index_j])))*sqrt(1.-alpha*alpha))    ;
      
      for(unsigned int ix=0 ;ix<3 ;ix++){
	domega1[j][ix] = - domega * dvac * (v[i][ix] - dotprod(v[i],v[index_j]) * v[index_j][ix]/norm2(v[index_j]));   // DER ANGULAR PART 1
	domega2[j][ix] =   domega * dvac * (v[index_j][ix] - dotprod(v[i],v[index_j]) * v[i][ix]/norm2(v[i]));          // DER ANGULAR PART 2
      }

      // TOTAL OF THE ANGULAR PART
      angtot += omega[j]*f[j];  // angular part multiplied by a weighting factor signifying radial contribution

    }
    

 // density swithing function   
    double drho=0.;                                                           
    double rho;                                                               // DENSITY SWITCHING FUNCTION
    rho = coordswitch.calculate(n,drho);					      // AND DERIVATIVES           


// Calculation of the CV = [(rho/n)*(f_ij*angular)] 
    Vector position = getPosition(i);
    Vector distrb = delta(Vector(0.,0.,0.),getPosition(i)); // calculate the distance of the atoms from the origin, center of the box
    double drb=distrb.modulo();   // scalar distance
    kernel(drb);        // calculate value of kernel function kval and its derivative dkval
  //double zpos = position[2];
  //kernel(zpos);   // calculate value of kernel function kval and its derivative dkval

    if(n>0.){           // THESE PARTS ARE USEFUL IN THE CALCULATION OF THE DERIVATIVES
      // NEW OPERATIONS TO SPEED UP
      S1 = drho * angtot;    
      S2 = -angtot/n;      
      S3 = rho*angtot/n;
      cv_val += S3*kval[2];     // SUM THE TOTAL CV 
    }else{
      S1 = 0.;
      S2 = 0.;
      S3 = 0.;
    }

   
    for( int j=0;j<GsmacShelllist.nn[i]; ++j) {                   // SUM OVER NEIGHBORS
      int index_j,start_j,end_j;
      int stride = (int)mols;
      index_j=GsmacShelllist.ni[i][j];                                      // TRUE INDEX OF THE J NEIGHBOR
      start_i= i + stride + stride;
      end_i= i + stride + stride + stride;
      start_j= index_j + stride;
      end_j=index_j + stride + stride;
      
      
      unsigned int m=3;
      for(unsigned int ix=0 ;ix<m ;ix++) {    ///  SMAC DERIVATIVES
	// CENTER
	double temp = ( (omega[j]  +  S2) * rho / n + S1) * df[j][ix] * kval[ix] + S3 * dkval[ix];
      	deriv[i][ix]          +=   temp;
      	deriv[index_j][ix]    -=   temp;
        // BOX DERIVATIVES
        dist=pbcDistance(getPosition(i),getPosition(index_j));    // DISTANCE BETWEEN THEM
 	for(unsigned jx=0;jx<m;jx++) {
      		 virial[ix][jx] += temp * dist[jx];
	}
	// START
	double temp1 = rho * f[j] * domega1[j][ix] / n * kval[ix];
	double temp2 = rho * f[j] * domega2[j][ix] / n * kval[ix];
	deriv[start_i][ix]    -=     temp2;     // ONLY ANGULAR PART
	deriv[start_j][ix]    +=     temp1;     // ONLY ANGULAR PART
	// END
	deriv[end_i][ix ]     +=     temp2;     // ONLY ANGULAR PART
	deriv[end_j][ix]      -=     temp1;     // ONLY ANGULAR PART

        // BOX DERIVATIVES
 	for(unsigned jx=0;jx<m;jx++) {
		 virial[ix][jx] -= temp2 * v[i][jx];
		 virial[ix][jx] += temp1 * v[index_j-mols][jx];
	}
      }
      
    }
  
    
  }
 
  comm.Sum(deriv);
  comm.Sum(cv_val);  
  comm.Sum(virial);  

 for(unsigned i=0;i<atoms;i++) {
   setAtomsDerivatives(i,deriv[i]);
 }
  
  setBoxDerivatives(virial);
  setValue(cv_val);
  
}
  
  
  
// Update and check list  
void GsmacShell::GsmacShell_newlist(vector<AtomNumber> &list1, GsmacShelllist_s &GsmacShelllist)
{
  Vector rij, test;
  double mod_rij;
  
  
  for(unsigned int j=0; j<list1.size(); ++j){
    GsmacShelllist.pos1[j] = getPosition(j);
    GsmacShelllist.pos2[j] = getPosition(j+list1.size());
  }
  
  unsigned int stride;
  unsigned int rank; 
  
  stride=comm.Get_size();  //Number of processes
  rank=comm.Get_rank(); //Rank of pr
  
    for(unsigned int j=rank; j<list1.size(); j+=stride) {     // sum over grid
      GsmacShelllist.nn[j]=0; //reset nlistsize
      for(unsigned int i=0; i<list1.size(); ++i) {                                           // sum over atoms
	if(i!=j){	  
	  rij = pbcDistance(GsmacShelllist.pos1[j],GsmacShelllist.pos2[i]);
	  mod_rij=rij.modulo2();

	  if (mod_rij < GsmacShelllist.rskin2){ //if distance < rskin
	    GsmacShelllist.ni[j][GsmacShelllist.nn[j]]= i + (unsigned int) list1.size(); //index in neighlist
	    GsmacShelllist.nn[j]++; //increment nn 
	  }
	}  
      }
    }
    
}

void GsmacShell::GsmacShell_checklist(vector<AtomNumber> &list1, GsmacShelllist_s &GsmacShelllist)
{
  unsigned int j; 
  double dr=(GsmacShelllist.rskin-GsmacShelllist.rcut)*0.5;
  Vector rij;
  
  for (j=0; j<list1.size(); ++j) { 
    //check position variations of center molecules
    rij = pbcDistance(getPosition(j),GsmacShelllist.pos1[j]); 
    if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
      GsmacShell_newlist(list1,GsmacShelllist); 
      break; 
    }
    //check position variations of neighbor molecules
    rij = pbcDistance(getPosition(j+list1.size()),GsmacShelllist.pos2[j]); 
    if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
      GsmacShell_newlist(list1,GsmacShelllist); 
      break; 
    }
  }
} 
    
double GsmacShell::norm2(Vector vect){           /// CALCULATE THE NORM OF A VECTOR
  return (vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]);
}

double GsmacShell::dotprod(Vector vect1,Vector vect2){           /// CALCULATE THE SCALAR PRODUCT OF TWO VECTORS
  return(vect1[0]*vect2[0]+vect1[1]*vect2[1]+vect1[2]*vect2[2]);
}
    
GsmacShell::~GsmacShell(){
}

  }
}

