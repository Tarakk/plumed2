/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_vesselbase_WeightedSumVessel_h
#define __PLUMED_vesselbase_WeightedSumVessel_h

#include <string>
#include <cstring>
#include <vector>
#include "VesselAccumulator.h"

namespace PLMD {
namespace vesselbase{

class WeightedSumVessel : public VesselAccumulator {
private:
  bool donorm;
protected:
/// We are normalizing the values
  void useNorm();
/// Add Some derivative to the weight 
  void addWeightDerivative( const unsigned&, const double& ); 
public:
  WeightedSumVessel( const VesselOptions& );
/// This retrieves data from action and calculates the average
  bool calculate( const unsigned& , const double& );
/// This does the final step of the calculation
  void finish( const double& tolerance );
/// This gets the weight
  virtual double getWeight( const unsigned&, bool& )=0;  
/// This adds the derivatives of the weight
  virtual void addDerivativesOfWeight( const unsigned& ){ plumed_merror("weights have no derivatives in this vessel"); }
/// This gets each value
  virtual double compute( const unsigned& , const unsigned& , const double& val, double& df )=0;
};

inline
void WeightedSumVessel::addWeightDerivative( const unsigned& ider, const double& der ){
  plumed_dbg_assert( ider<getNumberOfDerivatives(0) );
  addToBufferElement( ider+1, der );
}

}
}
#endif
