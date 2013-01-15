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
#ifndef __PLUMED_core_ActionAtomistic_h
#define __PLUMED_core_ActionAtomistic_h

#include "Action.h"
#include "tools/Tensor.h"
#include <vector>
#include <set>

namespace PLMD {

class Atoms;
class Pbc;
class PDB;

/// \ingroup MULTIINHERIT
/// Action used to create objects that access the positions of the atoms from the MD code
class ActionAtomistic :
  virtual public Action
  {

  std::vector<AtomNumber> indexes;         // the set of needed atoms
  std::set<AtomNumber>  unique;
  std::vector<Vector>   positions;       // positions of the needed atoms
  double                energy;
  Tensor                box;
  Pbc&                  pbc;
  Tensor                virial;
  std::vector<double>   masses;
  bool                  chargesWereSet;
  std::vector<double>   charges;

  std::vector<Vector>   forces;          // forces on the needed atoms
  double                forceOnEnergy;

  bool                  lockRequestAtoms; // forbid changes to request atoms

protected:
  Atoms&                atoms;

public:
/// Request an array of atoms.
/// This method is used to ask for a list of atoms. Atoms
/// should be asked for by number. If this routine is called
/// during the simulation, atoms will be available at the next step
/// MAYBE WE HAVE TO FIND SOMETHING MORE CLEAR FOR DYNAMIC
/// LISTS OF ATOMS
  void requestAtoms(const std::vector<AtomNumber> & a);
/// Get position of i-th atom
  const Vector & getPosition(int)const;
/// Get position of i-th atom
  const Tensor & getBox()const;
/// Get the array of all positions
  const std::vector<Vector> & getPositions()const;
/// Get energy
  const double & getEnergy()const;
/// Get mass of i-th atom
  double getMass(int i)const;
/// Get charge of i-th atom
  double getCharge(int i)const;
/// Get a reference to forces array
  std::vector<Vector> & modifyForces();
/// Get a reference to virial array
  Tensor & modifyVirial();
/// Get a reference to force on energy
  double & modifyForceOnEnergy();
/// Get number of available atoms
  unsigned getNumberOfAtoms()const{return indexes.size();};
/// Compute the pbc distance between two positions
  Vector pbcDistance(const Vector&,const Vector&)const;
/// Get the absolute index of an atom
  AtomNumber getAbsoluteIndex(int i)const;
/// Parse a list of atoms without a numbered keyword
  void parseAtomList(const std::string&key,std::vector<AtomNumber> &t);
/// Parse an list of atom with a numbred keyword
  void parseAtomList(const std::string&key,const int num, std::vector<AtomNumber> &t);
/// Get reference to Pbc
  const Pbc & getPbc() const;
public:

// virtual functions:

  ActionAtomistic(const ActionOptions&ao);
  ~ActionAtomistic();

  static void registerKeywords( Keywords& keys );

  void clearOutputForces();

/// N.B. only pass an ActionWithValue to this routine if you know exactly what you 
/// are doing.  The default will be correct for the vast majority of cases
  void   calculateNumericalDerivatives( ActionWithValue* a=NULL );

  void retrieveAtoms();
  void applyForces();
  void lockRequests();
  void unlockRequests();
  const std::set<AtomNumber> & getUnique()const;
/// Read in an input file containing atom positions and calculate the action for the atomic 
/// configuration therin
  void readAtomsFromPDB( const PDB& pdb );
};

inline
const Vector & ActionAtomistic::getPosition(int i)const{
  return positions[i];
}

inline
double ActionAtomistic::getMass(int i)const{
  return masses[i];
}

inline
double ActionAtomistic::getCharge(int i) const {
  if( !chargesWereSet ) error("charges were not passed to plumed");
  return charges[i];
}

inline
AtomNumber ActionAtomistic::getAbsoluteIndex(int i)const{
  return indexes[i];
}

inline
const std::vector<Vector> & ActionAtomistic::getPositions()const{
  return positions;
}

inline
const double & ActionAtomistic::getEnergy()const{
  return energy;
}

inline
const Tensor & ActionAtomistic::getBox()const{
  return box;
}

inline
std::vector<Vector> & ActionAtomistic::modifyForces(){
  return forces;
}

inline
Tensor & ActionAtomistic::modifyVirial(){
  return virial;
}

inline
void ActionAtomistic::clearOutputForces(){
  for(unsigned i=0;i<forces.size();++i)forces[i].zero();
  forceOnEnergy=0.0;
}


inline
double & ActionAtomistic::modifyForceOnEnergy(){
  return forceOnEnergy;
}

inline
const Pbc & ActionAtomistic::getPbc() const{
 return pbc;
}

inline
void ActionAtomistic::lockRequests(){
  lockRequestAtoms=true;
}

inline
void ActionAtomistic::unlockRequests(){
  lockRequestAtoms=false;
}

inline
const std::set<AtomNumber> & ActionAtomistic::getUnique()const{
  return unique;
}


}

#endif
