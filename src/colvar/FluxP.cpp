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

    //+PLUMEDOC COLVAR FLUXP 
    /*
      Calculates the number of molecules passing through an imaginary plane placed at a certain
      distance away from the interface. Slope of the #molecules vs time plot provides the flux.
      Use:
      # Define groups for FluxP
      GROUP ATOMS=1-5432:8 LABEL=urea        #all urea + water atoms
      GROUP ATOMS=5433-10703:3 LABEL=solv   #all water atoms
      FLUXP GROUPA=urea GROUPB=solv NINT=10 NZ=188 DISTP=1.0 STRIDE=1000 LABEL=fluxP
    */
    //+ENDPLUMEDOC
   
    class FluxP : public Colvar {
      bool isnotscaled;
      int  nbin;
      double  nint, distp;
      vector<unsigned> atomsRCl; //atoms in the Right (R) side of the crystal (C), left (l) of the plane
      vector<unsigned> atomsRCr; //atoms in the Right (R) side of the crystal (C), right (r) of the plane

      vector<unsigned> atomsLCl; //atoms in the Left  (L) side of the crystal (C), left (l) of the plane
      vector<unsigned> atomsLCr; //atoms in the Left  (L) side of the crystal (C), left (r) of the plane
      vector<AtomNumber> ga_lista,gb_lista,full_lista; // a: urea atoms, b: water atoms
      bool pbc;
      bool isFirstStep;
      int  storeHalfBin;
      int  stride;

      
    public:
      FluxP(const ActionOptions&);
      virtual void calculate();
      static void registerKeywords(Keywords& keys);
      ofstream fdbg;
    };

    PLUMED_REGISTER_ACTION(FluxP,"FLUXP")

    void FluxP::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUPA","Urea atoms");
      keys.add("atoms","GROUPB","water atoms)");
      keys.add("optional","NZ","Interface localization: zbin");
      keys.add("optional","NINT","Interface localization: density of solvent molecules");
      keys.add("compulsory","DISTP","distance from the right side of the plane");
    }

    FluxP::FluxP(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao),
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

     
      parse("DISTP",distp); //distance from the interface
 
      
      nint=0.0; //default values
      nbin=100;
      parse("NINT",nint); //interface boundary concentration
      parse("NZ",nbin); //z histogram bins
      log.flush(); //DBG
      checkRead();
      
      //request all atoms by combining two lists
      full_lista.insert( full_lista.end(), ga_lista.begin(), ga_lista.end() );
      full_lista.insert( full_lista.end(), gb_lista.begin(), gb_lista.end() );
      log.printf(" \n ");
      for(unsigned int i=0;i<full_lista.size();++i){
        if ( (i+1) % 25 == 0 ) log.printf("  \n");
        log.printf("  %d", full_lista[i].serial());
       }
      requestAtoms(full_lista);
      log.flush();       
      isFirstStep=true;
    // Add output component
      addComponent("fluxPL"); componentIsNotPeriodic("fluxPL");      
      addComponent("fluxPR"); componentIsNotPeriodic("fluxPR");      
      addComponent("zleft"); componentIsNotPeriodic("zleft");      
      addComponent("zright"); componentIsNotPeriodic("zright");      
    }

    
    // calculator
    void FluxP::calculate()    
    {
      double fluxPL,fluxP1L,fluxP2L;
      double fluxPR,fluxP1R,fluxP2R;
      vector<Vector> solv_x(gb_lista.size());      
      
      //Box size
      double LBC[3];
      for(int i=0;i<3;++i) LBC[i]=getBox()[i][i]; 
      

      //Histogram settings (for interface localization).............................................
      //histz-array allocation
      vector<int> histz(nbin,0.0);
      int nz=0;
      //bins width (time dependent)
      double dz=LBC[2]/nbin;
      double Vbin=LBC[0]*LBC[1]*dz; //Bin volume [nm^3]  

      for(int i=0; i<gb_lista.size(); i+=1){
          solv_x[i] = pbcDistance(Vector(0.,0.,0.),getPosition(ga_lista.size()+i));

        //impose PBC (useless on x and y for now!!!)
        if(solv_x[i][2]<0) solv_x[i][2]=solv_x[i][2]+LBC[2];
        nz=(int)(solv_x[i][2]/dz); //fill histogram
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

      //interface finder..............................................................................
      //Get the liquid-crystal interfaces in both sides
      double halfbin, ileft, iright, zleft, zright;
      halfbin=(int)(LBC[2]/(2*dz));
      int p=0;
      int pmone=0;

    if (isFirstStep) {
      isFirstStep=false;

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
      if(ileft<0) ileft=ileft+nbin; //pbc on left
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
      storeHalfBin=(ileft+iright)/2;
      while(ileft> nbin/2) {
         ileft=ileft-nbin; // effect of PBC in finding the interface
         storeHalfBin=(ileft+iright)/2;
      }

      if(ileft<0) ileft=ileft+nbin; //pbc on left
      zleft=dz*(ileft+1); //left interface coordinate
      zright=dz*(iright); //right interface coordinate
    
      // Calculation of fluxP...
      // calculate fluxP right side of the interface, checking the depletion.
      double ipr=zright+distp;
      if(ipr<0.0) ipr=ipr-LBC[2];

      // calculate number of molecules crossing the plane
      fluxP1R = 0.;
      for(unsigned int j=0;j<atomsRCl.size();j++){
        Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(atomsRCl[j]));      // position with resepect to origin
        if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];        //add pbc here
        if(pos[2]>=LBC[2]) pos[2]=pos[2]-LBC[2];     //add pbc here

        if(pbcDistance(Vector(0.,0.,ipr),pos)[2] > 0.){    // molecules that cross plane from left to right
          fluxP1R += 1.;
        }
      }
      atomsRCl.clear();

      // keep a track of molecules coming into the region defined
      fluxP2R = 0.;
      for(unsigned int j=0;j<atomsRCr.size();j++){
        Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(atomsRCr[j]));      // position with resepect to origin
        if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];          //add pbc here
        if(pos[2]>=LBC[2]) pos[2]=pos[2]-LBC[2];       //add pbc here

        if(pbcDistance(Vector(0.,0.,ipr),pos)[2] < 0.){  // molecules that cross plane from right to left
          fluxP2R += 1.;
        }
      }
      atomsRCr.clear();


      // Calculate fluxP in the left side of the crystal, checking growth rate.
      double ipl=zleft-distp;
      if(ipl<0.0) ipl=ipl+LBC[2];

      fluxP1L=0.;
      for(unsigned int j=0;j<atomsLCr.size();j++){
        Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(atomsLCr[j]));  // position with resepect to origin
        if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];      //add pbc here
        if(pos[2]>=LBC[2]) pos[2]=pos[2]-LBC[2];   //add pbc here

        if(pbcDistance(Vector(0.,0.,ipl),pos)[2]<0.){    // molecules that cross plane and move toward left
          fluxP1L += 1.;
        //log.printf("atoms in the left=%f %f %f\n",getAbsoluteIndex(j));
        }
      }
      atomsLCr.clear();

      // keep a track of molecules coming into the region defined
      fluxP2L = 0.;
      for(unsigned int j=0;j<atomsLCl.size();j++){
        Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(atomsLCl[j]));      // position with resepect to origin
        if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];        //add pbc here
        if(pos[2]>=LBC[2]) pos[2]=pos[2]-LBC[2];     //add pbc here

        if(pbcDistance(Vector(0.,0.,ipl),pos)[2] > 0.){          // molecules that cross plane from left to right
          fluxP2L += 1.;
        }
      }
      atomsLCl.clear();
    

      // Store atoms in left and right side of the plane
      for(unsigned int i=0;i<ga_lista.size();i+=1) {
         Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(i));      // position with resepect to origin
         if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];          //add pbc here
         if(pos[2]>=LBC[2]) pos[2]=pos[2]-LBC[2];       //add pbc here


       //Left side of the crystal
         if(pbcDistance(Vector(0.,0.,zleft),pos)[2]<0. && pbcDistance(Vector(0.,0.,ipl),pos)[2] >0.){          // molecules present within zleft <--> zleft-distp
             //log.printf("positions of right atoms %f \n", pos[2]);
             atomsLCr.push_back(i);             // storing atom indices into the vector atomsLCr
         }else if(pbcDistance(Vector(0.,0.,ipl),pos)[2]<0. && pbcDistance(Vector(0.,0.,(ipl-distp)),pos)[2] > 0.){  // molecules present within zleft-distp <--> zleft-2*distp
             //log.printf("positions of left atoms %f \n", pos[2]);
             atomsLCl.push_back(i);             //storing atom indices into the vector atomsLCl

       //Right side of the crystal
         }else if(pbcDistance(Vector(0.,0.,zright),pos)[2] > 0. && pbcDistance(Vector(0.,0.,ipr),pos)[2] < 0.){  // molecules present within zright <--> zright+distp
             atomsRCl.push_back(i);             // storing atom indices into the vector atomsRCl
         }else if(pbcDistance(Vector(0.,0.,ipr),pos)[2] > 0. && pbcDistance(Vector(0.,0.,(ipr+distp)),pos)[2] < 0.){ //molecules present within zright+distp <--> zright+2*distp 
             atomsRCr.push_back(i);             //storing atom indices into the vector atomsRCr
         }
      }

      // fluxP in the right side 
      fluxPR=fluxP1R-fluxP2R;

      // fluxP in the left side 
      fluxPL=fluxP2L-fluxP1L;
      

      getPntrToComponent("fluxPR")->set(fluxPR); //print out right side flux
      getPntrToComponent("fluxPL")->set(fluxPL); //print out left side flux

      getPntrToComponent("zleft")->set(zleft); //print out the interface coordinate
      getPntrToComponent("zright")->set(zright); //print out the interface coordinate
    }
  }  
}    
