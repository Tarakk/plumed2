/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   See http://www.plumed-code.org for more information.

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

   CmuMD method file, 
   see (and cite) Perego, Salvalaglio, Parrinello J. Chem. Phys. 142 144113 (2015) 
   http://scitation.aip.org/content/aip/journal/jcp/142/14/10.1063/1.4917200
   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
 
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

namespace PLMD{
  namespace colvar{

    //+PLUMEDOC COLVAR SHIFT 
    /*
    calculate the distance between the center of the simulation box and 
    center of the crystal; shift the entire simulation box to keep the crystal at the middle
    Use:
    # Define groups for shift
    GROUP ATOMS=1-10703 LABEL=all   #all water atoms
    SHIFT GROUP=all NSV=3 SOLUTE=5432 NST=8 NZ=188 NINT=5 STRIDE=1000 LABEL=shift
    */
    //+ENDPLUMEDOC
   
    class Shift : public Colvar {
      bool isnotscaled;
      bool isFirstStep;
      int storeHalfBin;
      int  N_st, N_sv, Na_sv_permol, Na_st_permol, Na_sv, Na_st, N_mol, com_sv, com_st, nbin;
      int Natot;
      double nint;
      
    public:
      Shift(const ActionOptions&);
      virtual void calculate();
      static void registerKeywords(Keywords& keys);
      double shift;
      int stride;
    };

    PLUMED_REGISTER_ACTION(Shift,"SHIFT")

    void Shift::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms involved in the calculation");
      keys.add("compulsory","NSV","Solvent atoms");
      keys.add("optional","SOLUTE","Solute tot atoms");
      keys.add("optional","NST","Solute atoms");
      keys.add("optional","NZ","Interface localization: zbin");
      keys.add("optional","NINT","Interface localization: int density");
      keys.add("optional","COMST","solute COM");
      keys.add("optional","COMSV","solvent COM");
      keys.add("compulsory","STRIDE","frequency at which shift will be done");
      keys.remove("NOPBC");
    }

    Shift::Shift(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao),
      //init bool parameters
      isnotscaled(false)
    {
      
      //Read atom group
      vector<AtomNumber> at_list;
      
      parseAtomList("GROUP",at_list);
      
      Na_sv_permol=1; //default
      parse("NSV",Na_sv_permol); //get number of atoms per molecule
      
      N_st=0; //default
      Na_st=0;
      Na_st_permol=1;
      
      parse("SOLUTE",Na_st); //get number of solute atoms
      parse("NST",Na_st_permol);
      
      //Solution numbers
      N_st=(int)(Na_st/Na_st_permol); //Number of solute atoms
      Na_sv=at_list.size()-Na_st; //Number of solvent atoms
      N_sv=(int)(Na_sv/Na_sv_permol); //Number of solvent molecules
      N_mol=N_sv+N_st; //Number of total molecules
      
      log.printf("Number of atoms:\tw %d\tu %d\n",Na_sv,Na_st);
      log.printf("Number of molecules:\ttot %d\t w %d\tu %d\n",N_mol,N_sv,N_st);
      log.flush();
      
      //COM flags
      com_sv=-1;
      com_st=-1;
      parse("COMSV",com_sv); 
      parse("COMST",com_st); 
      
      parse("STRIDE",stride); 
      
      nint=0.0; //default values
      nbin=100;
      parse("NINT",nint); //interface boundary concentration
      parse("NZ",nbin); //z histogram bins
      log.flush(); //DBG
      checkRead();
      addValueWithDerivatives(); 
      setNotPeriodic();
      
      requestAtoms(at_list);
      log.printf("  \n");
      log.flush();       
      isFirstStep=true;
    }

    
    // calculator
    void Shift::calculate()    
    {
      vector<Vector> com_solv(N_sv); 
      Vector diff;
      //Solvent position matrix allocation 
      vector<Vector> solve_x(Na_sv_permol);
      //Solvent mass array allocation 
      vector<double> solve_m(Na_sv_permol);
      
      //Solvent masses and total mass
      double M_sv=0.0;
      for(int i=0;i<Na_sv_permol;++i){
	solve_m[i]=getMass(Na_st+i); //the first Na_st are skipped
	M_sv += solve_m[i];
      }

      //Box size
      double LBC[3];
      for(int i=0;i<3;++i) LBC[i]=getBox()[i][i]; 
    
      
     // Restrain the crystal in the middle of the box
      Natot=Na_st+Na_sv;
      if(getStep()%stride==0) {
        // shift the coordinates of all atoms 
        for(int i=0; i<Natot; ++i){
           AtomNumber index;
           index = getAbsoluteIndex(i);
           Vector & ato (modifyPosition(index));
           ato += Vector(0.,0.,shift);
       
           // Wrap the coordinates, PBC on the right side
            if(ato[2]>=LBC[2]) ato[2]=ato[2]-LBC[2];
        }
      }  

      //Histogram settings (for interface localization)
      //histz-array allocation
      vector<int> histz(nbin,0.0);
      int nz=0;
      //bins width (time dependent)
      double dz=LBC[2]/nbin;
      double Vbin=LBC[0]*LBC[1]*dz; //Bin volume [nm^3]  
    
      //center of mass vector
      for(int i=0; i<N_sv; i+=1){
	com_solv[i].zero();
	//center of mass
	if (com_sv<0) {
	  solve_x[0] = getPosition(Na_st+i*Na_sv_permol);
	  for(int j=1;j<Na_sv_permol;++j){
	    solve_x[j] = getPosition(Na_st+i*Na_sv_permol+j);
	    diff = pbcDistance(solve_x[0],solve_x[j]);
	    com_solv[i] += solve_m[j]*diff;
	  }      
	  com_solv[i] = com_solv[i] / M_sv + solve_x[0]; 
	  //impose PBC on com (useless on x and y for now!!!)
	  if(com_solv[i][2]<0) com_solv[i][2]=com_solv[i][2]+LBC[2];
	  if(com_solv[i][2]>=LBC[2]) com_solv[i][2]=com_solv[i][2]-LBC[2];
	}else{
	  //no com
	  com_solv[i]=getPosition(Na_st+i*Na_sv_permol+com_sv);
	}
          nz=(int)(com_solv[i][2]/dz); //fill histogram
          histz[nz]+=1;
      }
 
      // Smooth histogram
      vector<double> smooth_histz(nbin,0.0);
       for(int i=0; i<nbin; i+=1){
         if (i==0) {
           smooth_histz[i]=(histz[i]+histz[i+1])/2.;
         } else if (i==(nbin-1)) {
           smooth_histz[i]=(histz[i]+histz[i-1])/2.;
         } else {
           smooth_histz[i]=(histz[i]+histz[i-1]+histz[i+1])/3.;
         }
         //log.printf("Histogram= %d %f\n",i,smooth_histz[i]); //useful to plot the solute distribution across the box
       }

       for(int i=0; i<nbin; i+=1){
            histz[i] = smooth_histz[i];
       }

      //communicate
      comm.Sum(histz);
      comm.Sum(com_solv);


     //Get the liquid-crystal interfaces
     //double halfbin, ileft, iright, zleft, zright;
     double halfbin, ileft, iright;
     double zleft, zright;

     //interface finder
      if (isFirstStep) {
          isFirstStep=false;
          halfbin=(int)(LBC[2]/(2*dz));
          int p=0;
          int pmone=0;
      

        //find the crystal if it's not at the half, it finds the crystal before halfbin exceeds the limits 
        //3 adjacent bins with water concentration < than nint/3
        while((histz[halfbin]+histz[halfbin+1]+histz[halfbin-1]) > nint*Vbin){
             p++;
             pmone=2*(p%2)-1;
             halfbin=halfbin+p*pmone; //Move through the bins
        }
     } else {
       halfbin=storeHalfBin;
     }
	
    //put halfbin inside the crystal volume (3 bins, WARNING!!! parameter dependent)
       ileft=halfbin;
	while(histz[ileft] < nint*Vbin){
	  ileft=ileft-1;
	  if(ileft<0) ileft=ileft+nbin; //pbc on left
	}
	
        iright=halfbin; //WARNING!!! parameter dependent
	if(iright>=nbin) iright=iright-nbin; //pbc on right
	while(histz[iright]< nint*Vbin){
	  iright=iright+1;
	  if(iright>=nbin) iright=iright-nbin; //pbc on right
	}
      
      // effect of PBC when the left interface is close to the left side
      // of the box, and eventually appears in the other side of the box due to PBC
      while(ileft> nbin/2) {
           ileft=ileft-nbin;
           storeHalfBin=(ileft+iright)/2;
      }
        storeHalfBin=(ileft+iright)/2;
	zleft=dz*(ileft+1); //left interface coordinate
	zright=dz*(iright); //right interface coordinate
      
       //calculate the shift
        double ccrystal=(zleft+zright)/2.0;
        shift = (LBC[2]/2.0) - ccrystal;   //calculation of the shift

     setValue(shift);     
   }  
  }    
}    
