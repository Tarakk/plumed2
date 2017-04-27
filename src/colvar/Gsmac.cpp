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
   
    class Gsmac : public Colvar{
      SwitchingFunction distanceswitch;      // NEED TO DEFINE MULTIPLE SWF
      SwitchingFunction coordswitch;      // NEED TO DEFINE MULTIPLE SWF
      
    private:


    public:

      vector<KernelFunctions> kernels;     //

      vector<AtomNumber> start;             // LIST OF MOLECULES NEEDED BY SMACK -> START
      vector<AtomNumber> end;               // LIST OF MOLECULES NEEDED BY SMACK -> END
      vector<AtomNumber> center1;            // LIST OF MOLECULES NEEDED BY SMACK -> CENTER1
      vector<AtomNumber> center2;            // LIST OF MOLECULES NEEDED BY SMACK -> CENTER2
      vector<AtomNumber> all_atoms;            // LIST OF MOLECULES NEEDED BY SMACK -> CENTER
      vector<double> angles;
      vector<double> width;

      unsigned int mols;                           //  total number of molecules
      unsigned int atoms;

      ofstream fdbg;
      ofstream rdbg;

      
      
      struct  Gsmaclist_s{       // STRUCTURE NEEDED TO CONSTRUCT THE VERLET LIST
	
	double rcut;             // CUT OFF OF THE LIST
	double rskin;            // SKIN TO CHECK IF THE LIST HAVE TO BE RECALCULATED
	double rskin2;            // SKIN SQUARED
	int  step;
	vector<int> nn;          // NUMBERS OF NEIGHBORS
	vector<vector<int> > ni; // LIST FOR ATOM I 
	vector<Vector> pos1;      // POSITIONS CENTER MOLECULE
	vector<Vector> pos2;      // POSITIONS NEIGHBORING MOLECULES
      } Gsmaclist;
      
      
      Gsmac(const ActionOptions&);              //    CONSTRUCTOR
      ~Gsmac();                                 //    DESTRUCTOR
      // active methods:
      static void registerKeywords( Keywords& keys );    // KEYWORDS
      virtual void calculate();                          // CALCULATE CV 
      //      void Gsmac_newlist(vector<AtomNumber> &list, vector<Vector> &grid, Gsmaclist_s &Gsmaclist);     // NEW VERLIST
      //      void Gsmac_checklist(vector<AtomNumber> &list, vector<Vector> &grid, Gsmaclist_s &Gsmaclist);   // CHECK THE VERLIST

      void Gsmac_newlist(vector<AtomNumber>&, Gsmaclist_s &Gsmaclist);     // NEW VERLIST
      void Gsmac_checklist(vector<AtomNumber>&, Gsmaclist_s &Gsmaclist);   // CHECK THE V
      double dotprod(Vector,Vector);
      double norm2(Vector);
      
    };
    
    PLUMED_REGISTER_ACTION(Gsmac,"GSMAC")
    
    void Gsmac::registerKeywords( Keywords& keys ){

      //  OrientationSphere::registerKeywords(keys);
      
      // IN DEFINING THE KEYWORDS I INTENTIONALLY DEFINE EVERYTHING HERE, WITHOUT CALLING ANY OTHER 
      // CLASS/FUNCTION/OBJECT SIMPLY BECAUSE i PREFER TO HAVE EVERYTHING UNDER CONTROL.
      Colvar::registerKeywords(keys);
      keys.add("atoms","CENTER1","the labels of the atoms acting as center of the molecules");

      keys.add("atoms","CENTER2","the labels of the atoms acting as center of the neighboring molecules");

      keys.add("atoms","START","the labels of the atoms acting as start of the intramolecular vector");

      keys.add("atoms","END","the labels of the atoms acting as end  of the intramolecular vector");

      keys.add("compulsory","ANGLES"," Angles that have to be used in the SMAC calculation (need to be associated with width)");

      keys.add("compulsory","N_ANGLES"," Number of angles that have to be used in the SMAC calculation (need to be associated with width)");

      keys.add("compulsory","WIDTH"," Width of the Gaussian used to calculate the angles distribution (need to be associated with angles)");

      keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
	       "The following provides information on the \\ref switchingfunction that are available. " 
	   "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords."); 
      
      keys.add("optional","SWITCH_COORD","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
	       "The following provides information on the \\ref switchingfunction that are available. " 
	       "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords."); 

      keys.add("compulsory","R_CUT","Verlet lists, default r_cut = R_max");
      keys.add("compulsory","R_SKIN","Verlet lists, default r_skin = 1.25 R_max");

      keys.remove("NOPBC");


    }
    
    Gsmac::Gsmac(const ActionOptions&ao):
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
      parseAtomList("CENTER2",center2);
      parseAtomList("START",start);
      parseAtomList("END",end);

      mols=center1.size();
      Gsmaclist.pos1.resize(mols);
      Gsmaclist.pos2.resize(mols);
      Gsmaclist.nn.resize(mols);
      Gsmaclist.ni.resize(mols);
      for(int i=0; i < mols ; i++){
	Gsmaclist.ni[i].resize(mols);
      }
      
      
      if(mols==0) error("no molecules specified");
      
      for(unsigned int i=0 ; i < mols ; i++){
	all_atoms.push_back(center1[i]);
      }
      for(unsigned int i=0 ; i < mols ; i++){
	all_atoms.push_back(center2[i]);
      }
      for(unsigned int i=0 ; i < mols ; i++){
	all_atoms.push_back(start[i]);
      }
      for(unsigned int i=0 ; i < mols ; i++){
	all_atoms.push_back(end[i]);
      }
      
      atoms=all_atoms.size();
      
      addValueWithDerivatives();       
      setNotPeriodic();
      requestAtoms(all_atoms);

      double r_cut,r_skin;
      parse("R_CUT",r_cut);
      parse("R_SKIN",r_skin);

      Gsmaclist.rcut=r_cut;
      Gsmaclist.rskin=r_skin;
      Gsmaclist.rskin2=r_skin*r_skin;
      
      checkRead();
      
      log<<"  contacts are counted with cutoff "<< distanceswitch.description()<<"\n";
      log<<"  correct density is evaluated above with cutoff "<< coordswitch.description()<<"\n";
      
      Gsmaclist.step=0;


      //CANCELLARE A FINE DEBUGGING
      //fdbg.open("dbg.dat");
      unsigned int rank;       
      rank=comm.Get_rank(); //Rank of pr
      stringstream A;
      //A<< "rdbg.dat."<< (int)rank;
      //rdbg.open(A.str().c_str());
      //fdbg << "DEBUG FILE \n";  //DBG
      //fdbg.flush(); //DBG
      //log.printf("dbg.dat open  \n");

    }
    
    
void Gsmac::calculate()
{
  
  double cv_val;             // CV
  cv_val=0;
  
  Tensor virial;              // VIRIAL
  vector<Vector> deriv(getNumberOfAtoms());  // DERIVATIVES
  
  if(Gsmaclist.step==0) Gsmac_newlist(center1,Gsmaclist); // CHECK IF NEIGH LIST HAVE TO BE CONSTRUCTED
  Gsmaclist.step=1; Gsmac_checklist(center1,Gsmaclist);   // CALL NEIGH LIST IN CASE
  
  unsigned int stride=comm.Get_size();                        // SET THE PARALLELIZATION VARIABLES
  unsigned int rank=comm.Get_rank();
  //      if(serial){                                             // SERIAL
  //	stride=1;
  //	rank=0;
  //      }else{                                                  // PARALLEL
  stride=comm.Get_size();
  rank=comm.Get_rank();
  //      }

  //parallel DBG
  
  //open rank debug file DBG
  
  vector<Vector> v(getNumberOfAtoms());                          // 
  for(unsigned int i = rank; i < mols; i += stride) {                   // DEFINE THE INTRAMOLECULAR VECTOR
    int start=i + (int) mols + (int) mols;
    int end=i + (int) mols + (int) mols + (int) mols;
    v[i]=pbcDistance(getPosition(start),getPosition(end));        // THIS IS DONE ONCE PER STEPS FOR ALL MOLECULES
  }
  
  comm.Sum(v);  // MERGING UP

  for(unsigned int i=rank;i<mols;i+=stride) {                   // SUM OVER MOLECULES
    double n,angtot;
    int start_i,end_i;
    double S1,S2;
    start_i=start[i].serial();
    end_i=end[i].serial();

    
    vector<double> f(Gsmaclist.nn[i]);                            // SWITCHING FUNCTION
    vector<double> omega(Gsmaclist.nn[i]);                            // SWITCHING FUNCTION
    vector<Vector> domega1(Gsmaclist.nn[i]);                      // ANGULAR PART1
    vector<Vector> domega2(Gsmaclist.nn[i]);                      // ANGULAR PART2
    vector<Vector> df(Gsmaclist.nn[i]);                      // ANGULAR PART2
    
    n=0.;
    angtot=0.;
    
    Vector dist;

    for( int j=0; j < Gsmaclist.nn[i]; ++j) {                   // SUM OVER NEIGHBORS
      int index_j;
      index_j=Gsmaclist.ni[i][j]-mols;                                      // TRUE INDEX OF THE J NEIGHBOR
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

      angtot += omega[j]*f[j];                                                                                       // TOTAL OF THE ANGULAR PART

    }
    

    
    double drho=0.;                                                           
    double rho;                                                               // DENSITY SWITCHING FUNCTION
    rho = coordswitch.calculate(n,drho);					      // AND DERIVATIVES           
    
    if(n>0.){           // THESE PARTS ARE USEFUL IN THE CALCULATION OF THE DERIVATIVES
      // NEW OPERATIONS TO SPEED UP
      S1 = drho * angtot;    
      S2 = -angtot / n;      
      cv_val += rho*angtot/(n);     // SUM THE TOTAL CV 
    }else{
      S1 = 0.;
      S2 = 0.;
    }

   
    for( int j=0;j<Gsmaclist.nn[i]; ++j) {                   // SUM OVER NEIGHBORS
      int index_j,start_j,end_j;
      int stride = (int)mols;
      index_j=Gsmaclist.ni[i][j];                                      // TRUE INDEX OF THE J NEIGHBOR
      start_i= i + stride + stride;
      end_i= i + stride + stride + stride;
      start_j= index_j + stride;
      end_j=index_j + stride + stride;
      
      
      unsigned int m=3;
      for(unsigned int ix=0 ;ix<m ;ix++) {    ///  SMAC DERIVATIVES
	// CENTER
	double temp = ( (omega[j]  +  S2) * rho / n + S1) * df[j][ix];
      	deriv[i][ix]          +=   temp ;
      	deriv[index_j][ix]    -=   temp;
        // BOX DERIVATIVES
        dist=pbcDistance(getPosition(i),getPosition(index_j));    // DISTANCE BETWEEN THEM
 	for(unsigned jx=0;jx<m;jx++) {
      		 virial[ix][jx] += temp * dist[jx];
	}
	// START
	double temp1 = rho * f[j] * domega1[j][ix] / n;
	double temp2 = rho * f[j] * domega2[j][ix] / n;
	deriv[start_i][ix]    -=     temp2;       // ONLY ANGULAR PART
	deriv[start_j][ix]    +=     temp1;     // ONLY ANGULAR PART
	// END
	deriv[end_i][ix ]     +=     temp2;     // ONLY ANGULAR PART
	deriv[end_j][ix]      -=     temp1;     // ONLY ANGULAR PART

        // BOX DERIVATIVES
 	for(unsigned jx=0;jx<m;jx++) {
		 virial[ix][jx] -= temp2 * v[i][jx] ;
		 virial[ix][jx] += temp1 * v[index_j-mols][jx] ;
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
  
  
  
  
void Gsmac::Gsmac_newlist(vector<AtomNumber> &list1, Gsmaclist_s &Gsmaclist)
{
  Vector rij, test;
  double mod_rij;
  
  
  for(unsigned int j=0; j<list1.size(); ++j){
    //test=getPosition(j);
    Gsmaclist.pos1[j] = getPosition(j);
    Gsmaclist.pos2[j] = getPosition(j+list1.size());
  }
  
  unsigned int stride;
  unsigned int rank; 
  
  stride=comm.Get_size();  //Number of processes
  rank=comm.Get_rank(); //Rank of pr
  
    for(unsigned int j=rank; j<list1.size(); j+=stride) {     // sum over grid
      Gsmaclist.nn[j]=0; //reset nlistsize
      for(unsigned int i=0; i<list1.size(); ++i) {                                           // sum over atoms
	if(i!=j){	  
	  rij = pbcDistance(Gsmaclist.pos1[j],Gsmaclist.pos2[i]);
	  mod_rij=rij.modulo2();

	  if (mod_rij < Gsmaclist.rskin2){ //if distance < rskin
	    Gsmaclist.ni[j][Gsmaclist.nn[j]]= i + (unsigned int) list1.size(); //index in neighlist
	    Gsmaclist.nn[j]++; //increment nn 
	  }
	}  
      }
    }
    
    //    comm.Sum(Gsmaclist.nn);
    //    comm.Sum(Gsmaclist.ni[0]);
}

void Gsmac::Gsmac_checklist(vector<AtomNumber> &list1, Gsmaclist_s &Gsmaclist)
{
  unsigned int j; 
  double dr=(Gsmaclist.rskin-Gsmaclist.rcut)*0.5;
  Vector rij;
        
  //fdbg << dulist.step << endl;
  //fdbg.flush();
  
  
  for (j=0; j<list1.size(); ++j) { 
    //check position variations of center molecules
    rij = pbcDistance(getPosition(j),Gsmaclist.pos1[j]); 
    if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
      //	id++; //rebuild list counter
      //fdbg << "rebuild list"<< endl;
      Gsmac_newlist(list1,Gsmaclist); 
      break; 
    }
    //check position variations of neighbor molecules
    rij = pbcDistance(getPosition(j+list1.size()),Gsmaclist.pos2[j]); 
    if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) { 
      Gsmac_newlist(list1,Gsmaclist); 
      break; 
    }
  }
} 
    
double Gsmac::norm2(Vector vect){           /// CALCULATE THE NORM OF A VECTOR
  return (vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]);
}

double Gsmac::dotprod(Vector vect1,Vector vect2){           /// CALCULATE THE SCALAR PRODUCT OF TWO VECTORS
  return(vect1[0]*vect2[0]+vect1[1]*vect2[1]+vect1[2]*vect2[2]);
}
    
Gsmac::~Gsmac(){
}

  }
}

