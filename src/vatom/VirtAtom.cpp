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
#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
 
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

namespace PLMD{
  namespace vatom{

    //+PLUMEDOC VATOM VIRTATM 
    /*
      Calculates the position of the VIRTUAL ATOM at the interface 
      Use:
      # Define groups for VirtAtom
      GROUP ATOMS=1-5432:8 LABEL=urea        #all urea + water atoms
      GROUP ATOMS=5433-10703:3 LABEL=solv   #all water atoms
      VIRTATM GROUPA=urea GROUPB=solv NINT=10 NZ=188 LABEL=v1
    */
    //+ENDPLUMEDOC
   
    class VirtAtom : public ActionWithVirtualAtom {
      bool nopbc;
      bool isnotscaled;
      int  nbin;
      double  nint;
      vector<unsigned> atomsR; //atoms in the Right (R) side of the crystal (C), left (l) of the plane
      vector<unsigned> atomsL; //atoms in the Left  (L) side of the crystal (C), left (r) of the plane
      vector<AtomNumber> ga_lista,gb_lista; // a: urea atoms, b: water atoms
      bool           pbc;

      
    public:
      explicit VirtAtom(const ActionOptions&);
      virtual void calculate();
      static void registerKeywords(Keywords& keys);
      ofstream fdbg;
    };

    PLUMED_REGISTER_ACTION(VirtAtom,"VIRTATM")

    void VirtAtom::registerKeywords(Keywords& keys){
      ActionWithVirtualAtom::registerKeywords(keys);
      keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating Distances");
      keys.add("atoms","GROUPA","Urea atoms");
      keys.add("atoms","GROUPB","water atoms)");
      keys.add("optional","NZ","Interface localization: zbin");
      keys.add("optional","NINT","Interface localization: density of solvent molecules");
    }

    VirtAtom::VirtAtom(const ActionOptions&ao):
      Action(ao),
      ActionWithVirtualAtom(ao),
      isnotscaled(false),
      nopbc(true)
    {
      parseFlag("NOPBC",nopbc);

      
      parseAtomList("GROUPA",ga_lista); //atoms in urea molecules
      parseAtomList("GROUPB",gb_lista); //atoms in solvent (water) molecules
       
          
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
    }

    
    // calculator
    void VirtAtom::calculate()    
    {
      vector<Tensor> deriv(getNumberOfAtoms());
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
 
       double mass=1.0;
       // Z-coordinate of the virtual atom 
       setPosition   (Vector(LBC[0]/2,LBC[1]/2,zright));
       setMass(mass);
       setAtomsDerivatives(deriv);
    }
  }  
}    
