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
      distance away from the interface. 
      Use:
      # Define groups for FluxP
      GROUP ATOMS=1-5432:8 LABEL=urea        #all urea + water atoms
      GROUP ATOMS=5433-10703:3 LABEL=solv   #all water atoms
      FLUXP GROUPA=urea GROUPB=solv NINT=10 NZ=188 DISTP=0.05 LABEL=fluxP
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
      vector<AtomNumber> ga_lista,gb_lista; // a: urea atoms, b: water atoms
      bool           pbc;

      
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

     
      parse("DISTP",distp); //distance from the left side of the interface
 
      
      nint=0.0; //default values
      nbin=100;
      parse("NINT",nint); //interface boundary concentration
      parse("NZ",nbin); //z histogram bins
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
    //addComponent("intleft"); componentIsNotPeriodic("inleft");      
    //addComponent("intright"); componentIsNotPeriodic("intright");      
      addComponent("fluxPL"); componentIsNotPeriodic("fluxPL");      
      addComponent("fluxPR"); componentIsNotPeriodic("fluxPR");      
      addComponent("fluxP"); componentIsNotPeriodic("fluxP");      
    }

    
    // calculator
    void FluxP::calculate()    
    {
      double fluxP;
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
       log.printf("zleft=%f\n",zleft);
       log.printf("zright=%f\n",zright);
     //double pos_left=zleft;
     //double pos_right=zright;
 
     //getPntrToComponent("intleft")->set(pos_left); //print out the interface coordinate
     //getPntrToComponent("intright")->set(pos_right); //print out the interface coordinate
     

      // Calculation of fluxP
      // calculate fluxP right side of the interface, checking the depletion.................................................................
      // calculate number of molecules crossing the plane
      fluxP1R = 0.;
      for(unsigned int j=0;j<atomsRCl.size();j++){
        Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(atomsRCl[j]));      // position with resepect to origin
        if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];    //add pbc here
        //if(pos[2]>=LBC[2]) pos[2]=pos[2]-LBC[2];                           //add pbc here
        double distancez = pos[2];
        if((distancez >= zright+(distp)) ){      //imaginary plane has been placed distp distance away from the interface
          fluxP1R += 1.;
        //log.printf("atoms in the left=%f %f %f\n",getAbsoluteIndex(j));
        }
      }
      atomsRCl.clear();

      // keep a track of molecules coming into the region defined
      fluxP2R = 0.;
      for(unsigned int j=0;j<atomsRCr.size();j++){
        Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(atomsRCr[j]));      // position with resepect to origin
        if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];                          //add pbc here
        //if(pos[2]>=LBC[2]) pos[2]=pos[2]-LBC[2];                     //add pbc here
        double distancez = pos[2];
        if((distancez <= zright+(distp)) ){                      //imaginary plane has been placed distp distance away from the interface
          fluxP2R += 1.;
        }
      }
      atomsRCr.clear();


      // Calculate fluxP in the left side of the crystal (check the growth rate)...............................................................
      fluxP1L=0.;
      for(unsigned int j=0;j<atomsLCr.size();j++){
        Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(atomsLCr[j]));      // position with resepect to origin
        if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];    //add pbc here
        //if(pos[2]>=LBC[2]) pos[2]=pos[2]-LBC[2];                           //add pbc here
        double distancez = pos[2];
        if((distancez <= zleft-(distp)) ){      //imaginary plane has been placed distp distance away from the interface
          fluxP1L += 1.;
        //log.printf("atoms in the left=%f %f %f\n",getAbsoluteIndex(j));
        }
      }
      atomsLCr.clear();

      // keep a track of molecules coming into the region defined
      fluxP2L = 0.;
      for(unsigned int j=0;j<atomsLCl.size();j++){
        Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(atomsLCl[j]));      // position with resepect to origin
        if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];                          //add pbc here
        //if(pos[2]>=LBC[2]) pos[2]=pos[2]-LBC[2];                     //add pbc here
        double distancez = pos[2];
        if((distancez >= zleft-(distp)) ){                      //imaginary plane has been placed distp distance away from the interface
          fluxP2L += 1.;
        }
      }
      atomsLCl.clear();

      //Left side of the crystal
      // Store atoms in left and right side of the plane
      for(unsigned int i=0;i<ga_lista.size();i+=1) {
        Vector pos = pbcDistance(Vector(0.,0.,0.),getPosition(i));      // position with resepect to origin
        if(pos[2]<=0.0) pos[2]=pos[2]+LBC[2];                           //add pbc here
        //if(pos[2]>=LBC[2]) pos[2]=pos[2]-LBC[2];                      //add pbc here
        double distancez = pos[2];
        if((distancez <= zleft) && (distancez >= zleft-(distp))){          // molecules that are present within zleft <--> zleft+distp
            atomsLCr.push_back(i);             // storing atom indices into the vector atomsL
          //log.printf("atoms in the left=%d \n ",ga_lista[i].serial());
        }else if(distancez <= zleft-(distp) && (distancez >= zleft-2*(distp))){     // molecules that are present inside a block right after the plane
            atomsLCl.push_back(i);             //storing atom indices into the vector atomsR

      //Right side of the crystal
      //store atoms in left and right side of the plane 
          //log.printf("atoms in the right=%d \n ",ga_lista[i].serial());
        }else if((distancez >= zright) && (distancez <= zright+(distp))){          // molecules that are present within zleft <--> zleft+distp
            atomsRCl.push_back(i);             // storing atom indices into the vector atomsL
          //log.printf("atoms in the left=%d \n ",ga_lista[i].serial());
        }else if(distancez >= zright+(distp) && (distancez <= zright+2*(distp))){     // molecules that are present inside a block right after the plane
            atomsRCr.push_back(i);             //storing atom indices into the vector atomsR
          //log.printf("atoms in the right=%d \n ",ga_lista[i].serial());
        }
      }

      // Total fluxP in the right side fluxPR= fluxP1-fluxP2
      fluxPR=fluxP2R-fluxP1R;

      // Total fluxP in the left side fluxPL= fluxP1L-fluxP2L
      fluxPL=fluxP2L-fluxP1L;
      
      // difference in the fluxPes 
      fluxP=abs(abs(fluxPL)-abs(fluxPR));

      getPntrToComponent("fluxPR")->set(fluxPR); //print out the interface coordinate
      getPntrToComponent("fluxPL")->set(fluxPL); //print out the interface coordinate

      //setValue(fluxPR);
      getPntrToComponent("fluxP")->set(fluxP);
      //setBoxDerivativesNoPbc();
    }
  }  
}    
