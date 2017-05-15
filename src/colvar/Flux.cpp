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

    //+PLUMEDOC COLVAR FLUX 
    /*
      Calculates the number of molecules passing through an imaginary plane placed at a certain
      distance away from the interface. 
      Use:
      # Define groups for Flux
      GROUP ATOMS=1-5432:8 LABEL=urea        #all urea + water atoms
      GROUP ATOMS=5433-10703:3 LABEL=solv   #all water atoms
      FLUX GROUPA=urea GROUPB=solv NINT=10 NZ=188 zi=0.5 zf=1.0 LABEL=flux
    */
    //+ENDPLUMEDOC
   
    class Flux : public Colvar {
      bool isnotscaled;
      int  nbin;
      unsigned stride_;
      double  nint;
      double zi,zf,dzfi;
      vector<AtomNumber> ga_lista,gb_lista; // a: urea atoms, b: water atoms
      bool           pbc;

      
    public:
      Flux(const ActionOptions&);
      virtual void calculate();
      static void registerKeywords(Keywords& keys);
      ofstream fdbg;
    };

    PLUMED_REGISTER_ACTION(Flux,"FLUX")

    void Flux::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUPA","Urea atoms");
      keys.add("atoms","GROUPB","water atoms)");
      keys.add("optional","NZ","Interface localization: zbin");
      keys.add("optional","NINT","Interface localization: density of solvent molecules");
      keys.add("compulsory","zi","z initial for the slab");
      keys.add("compulsory","zf","z final for the slab, dzfi=|zf-zi|");
      keys.add("compulsory","STRIDE","the frequency with which the flux should be calculated.");
    }

    Flux::Flux(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao),
      //init bool parameters
      isnotscaled(false),
      pbc(true)
    {
      bool nopbc=!pbc;
      parseFlag("NOPBC",nopbc);
      pbc=!nopbc;

      
       parseAtomList("GROUPA",ga_lista); //atoms in urea molecules
       parseAtomList("GROUPB",gb_lista); //atoms in solvent (water) molecules
       
       log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
       log.printf("  first group:\n");
       for(unsigned int i=0;i<ga_lista.size();++i){
        if ( (i+1) % 25 == 0 ) log.printf("  \n");
        log.printf("  %d", ga_lista[i].serial());
       }
       log.printf("  \n  second group:\n");
       for(unsigned int i=0;i<gb_lista.size();++i){
        if ( (i+1) % 25 == 0 ) log.printf("  \n");
        log.printf("  %d", gb_lista[i].serial());
       }

     
      parse("zi",zi); //distance from the left side of the interface
      parse("zf",zf); //distance from the left side of the interface
   
      dzfi=abs(zf-zi); 
      
      nint=0.0; //default values
      nbin=100;
      parse("NINT",nint); //interface boundary concentration
      parse("NZ",nbin); //z histogram bins

      parse("STRIDE",stride_); //frequency of calculating flux
      if(stride_<=0) error("frequency for flux calculation is nonsensical");

      log.flush(); //DBG
      checkRead();
      
      //request all atoms by combining two lists
      vector<AtomNumber> full_lista;      
      full_lista.insert( full_lista.end(), ga_lista.begin(), ga_lista.end() );
      full_lista.insert( full_lista.end(), gb_lista.begin(), gb_lista.end() );
      log.printf(" \n ");
      for(unsigned int i=0;i<full_lista.size();++i){
        if ( (i+1) % 25 == 0 ) log.printf("  \n");
        log.printf("  %d", full_lista[i].serial());
       }
      requestAtoms(full_lista);
      log.flush();       

    // Add output component
      addComponent("zleft"); componentIsNotPeriodic("zleft");      
      addComponent("zright"); componentIsNotPeriodic("zright");      
      addComponent("fluxL"); componentIsNotPeriodic("fluxL");      
      addComponent("fluxR"); componentIsNotPeriodic("fluxR");      
    }

    
    // calculator
    void Flux::calculate()    
    {
      double flux;
      double fluxL;
      double fluxR;
      vector<Vector> solut_x(ga_lista.size());      
      vector<Vector> solv_x(gb_lista.size());      
      vector<Vector> v(ga_lista.size()); //define a vector that stores velocities of the particles     
      
      //Box size
      double LBC[3];
      for(int i=0;i<3;++i) LBC[i]=getBox()[i][i]; 

      // printing velocities of atoms
    //for(int i=0; i<ga_lista.size(); i+=1){
    //    v[i] = getVelocity(i);
    //    log.printf("velocities of atoms=%f %f %f\n",v[i][0],v[i][1],v[i][2]);
    //}
      

      //Histogram settings (for interface localization).............................................
      //histz-array allocation
      vector<int> histz(nbin,0.0);
      int nz=0;
      //bins width (time dependent)
      double dz=LBC[2]/nbin;
      double Vbin=LBC[0]*LBC[1]*dz; //Bin volume [nm^3]  

      // Flux histogram, solute distribution
      vector<int> histu(nbin,0.0);
      for(int i=0; i<ga_lista.size(); i+=1){
            solut_x[i] = pbcDistance(Vector(0.,0.,0.),getPosition(i));

          //impose PBC (useless on x and y for now!!!)
          if(solut_x[i][2]<0) solut_x[i][2]=solut_x[i][2]+LBC[2];
          nz=(int)(solut_x[i][2]/dz); //fill histogram
          histu[nz]+=1.;
          //histu[nz]+=getVelocity(i)[2];
      }

      // Smooth histogram
      vector<double> smooth_histu(nbin,0.0);
      for(int i=0; i<nbin; i+=1){
         if (i==0) {
           smooth_histu[i]=(histu[i]+histu[i+1])/2.;
         } else if (i==(nbin-1)) {
           smooth_histu[i]=(histu[i]+histu[i-1])/2.;
         } else {
           smooth_histu[i]=(histu[i]+histu[i-1]+histu[i+1])/3.;
         }
       //log.printf("Histogram= %d %f\n",i,smooth_histu[i]);
      }
      

      // Solvent distribution
      for(int i=0; i<gb_lista.size(); i+=1){
          solv_x[i] = pbcDistance(Vector(0.,0.,0.),getPosition(ga_lista.size()+i));

        //impose PBC (useless on x and y for now!!!)
        if(solv_x[i][2]<0) solv_x[i][2]=solv_x[i][2]+LBC[2];
        nz=(int)(solv_x[i][2]/dz); //fill histogram
        histz[nz]+=1;
      }

      //interface finder..............................................................................
      //Get the liquid-crystal interfaces in both sides
      double halfbin, ileft, iright, zleft, zright;
      halfbin=(int)(LBC[2]/(2*dz));
      int p=0;
      int pmone=0;

      while((histz[halfbin]+histz[halfbin+1]+histz[halfbin-1]) > nint*Vbin){
         p++;
         pmone=2*(p%2)-1;
         halfbin=halfbin+p*pmone; //Move through the bins
       }
       //put halfbin inside the crystal volume (3 bins, WARNING!!! parameter dependent)
       ileft=halfbin;
       while(histz[ileft] < nint*Vbin){
         ileft=ileft-1;
         if(ileft<0) ileft=ileft+nbin; //pbc on left
       }

       iright=ileft+10; //WARNING!!! parameter dependent
       if(iright>=nbin) iright=iright-nbin; //pbc on right
       while(histz[iright]< nint*Vbin){
         iright=iright+1;
         if(iright>=nbin) iright=iright-nbin; //pbc on right
       }

       zleft=dz*(ileft+1); //left interface coordinate
       zright=dz*(iright); //right interface coordinate
 
     
    fluxL = 0.;
    fluxR = 0.;

      // Calculation of flux
      for(unsigned int i=0;i<ga_lista.size();i+=1) {
        Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(i));      // position with resepect to origin
        if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];                           //add pbc here
        double distancez = pos[2];
        AtomNumber index;
	index = getAbsoluteIndex(i);
        //Left side of the crystal
        if((distancez <= zleft-zi) && (distancez >= zleft-(zi+zf))){          // molecules that are present within zleft <--> zleft+distp
          fluxL += getVelocity(index)[2];
        //Right side of the crystal
        }else if((distancez >= zright+zi) && (distancez <= zright+(zi+zf))){          // molecules that are present within zleft <--> zleft+distp
          fluxR += getVelocity(index)[2];
        }
      }


    if(getStep()%stride_==0) {

    // rescale the velocity of the particles present in the right side of the crystal
    double k=-fluxL/fluxR; //factor to rescale velocity
    log.printf("k=%f\n",k);
    for(unsigned int i=0;i<ga_lista.size();i+=1) {
      Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(i));      // position with resepect to origin
      if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];                           //add pbc here
        double distancez = pos[2];
        AtomNumber index;
	index = getAbsoluteIndex(i);
        //Right side of the crystal
        if((distancez >= zright+zi) && (distancez <= zright+(zi+zf))){          
           //log.printf("before - verocities right side=%f\n",getVelocity(index)[2]);
           rescaleVelocity(index,Vector(1.,1.,k));
        }
    } 

    fluxR = 0.;

      // Calculation of flux
      for(unsigned int i=0;i<ga_lista.size();i+=1) {
        Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(i));      // position with resepect to origin
        if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];                           //add pbc here
        double distancez = pos[2];
        AtomNumber index;
	index = getAbsoluteIndex(i);
        //Right side of the crystal
        if((distancez >= zright+zi) && (distancez <= zright+(zi+zf))){          // molecules that are present within zleft <--> zleft+distp
          //fluxR += getVelocity(getAbsoluteIndex(i))[2];
          fluxR += getVelocity(index)[2];
        }
      }

    }


      fluxL=fluxL/dzfi; //dividing the flux by the z-box distance
      fluxR=fluxR/dzfi;

      // difference in the fluxes 
      //flux=abs(abs(fluxL)-abs(fluxR));

      getPntrToComponent("zleft")->set(zleft); //print out the left interface coordinate
      getPntrToComponent("zright")->set(zright); //print out the right interface coordinate

      getPntrToComponent("fluxR")->set(fluxR); //print out the interface coordinate
      getPntrToComponent("fluxL")->set(fluxL); //print out the interface coordinate

      //setValue(flux);
      //setBoxDerivativesNoPbc();

    }
  }  
}    
