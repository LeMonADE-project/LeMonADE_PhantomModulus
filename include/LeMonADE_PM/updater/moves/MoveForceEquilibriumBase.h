/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
    ooo                        |
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#ifndef LEMONADE_PM_UPDATER_MOVES_MOVEFORCEEQUILIBRIUMBASE_H
#define LEMONADE_PM_UPDATER_MOVES_MOVEFORCEEQUILIBRIUMBASE_H

#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveForceEquilibriumBase
 *
 * @brief Base class for all force equilibration-moves. 
 *
 *
 * @tparam <SpecializedMove> name of the specialized move.
 *
 **/
/*****************************************************************************/
template <class SpecializedMove>
class MoveForceEquilibriumBase:public MoveBase
{
 public:
	//! Returns the index of the Vertex (monomer) in the graph which searches for a partner
	uint32_t getIndex() const {return index;}

	//! Returns the direction of the Vertex (monomer) in the graph which should be moved
	const VectorDouble3& getShiftVector() const { return ShiftVector;}
	
	
	//here come the functions that are implemented by the specialization
	template <class IngredientsType> void init(const IngredientsType& ingredients);
	template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
	template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, VectorDouble3 dir);

	template <class IngredientsType> void check(const IngredientsType& ingredients);
	template <class IngredientsType> void apply(IngredientsType& ingredients);

 protected:
	/**
	 * @brief Set the index of the Vertex (monomer) in the graph which should be moved
	 * @param i the index in the graph
	 */
	void setIndex(uint32_t i) {index=i;}

	/**
	 * @brief Set the move direction of the Vertex (monomer) which should be moved
	 * @param dir The displacement in direction on the Cartesian-space.
	 */
	void setShiftVector(const VectorDouble3& ShiftVector_) {ShiftVector=ShiftVector_;}

	/**
	 * @brief Set the move direction of the Vertex (monomer) which should be moved
	 *
	 * @param dx The displacement in x-direction in the Cartesian space
	 * @param dy The displacement in y-direction in the Cartesian space
	 * @param dz The displacement in z-direction in the Cartesian space
	 */
	void setShiftVector(const double dx, const double dy, const double dz)
	{
		ShiftVector.setAllCoordinates(dx,dy,dz);
	}
public:
    void setFilename(std::string filename_){static_cast<SpecializedMove*>(this)->setFilename(filename_);}
    //! get the filename for the force extension data 
    std::string const getFilename(){static_cast<SpecializedMove*>(this)->getFilename();}
    
    //! set the relaxation parameter for the cross link
    void setRelaxationParameter(double relaxationChain_){static_cast<SpecializedMove*>(this)->setRelaxationParameter(relaxationChain_);}
    //! get the relaxation parameter for the cross link 
    double getRelaxationParameter(){static_cast<SpecializedMove*>(this)->getRelaxationParameter();} 


	//! Random Number Generator (RNG)
	RandomNumberGenerators randomNumbers;


 private:

	//! Index of the Vertex (monomer) in the graph which searches for a partner
	uint32_t index;

	//! Direction for the move
	VectorDouble3 ShiftVector;
	
};

////////////////////////////////////////////////////////////////////////////////
// implementation of the members
////////////////////////////////////////////////////////////////////////////////
/*****************************************************************************/
/**
 * @brief Initialize the move. Done by the SpecializedMove.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter
 *
 * @tparam <SpecializedMove> name of the specialized move.
 **/
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType>
void MoveForceEquilibriumBase<SpecializedMove>::init(const IngredientsType& ingredients)
{
  static_cast<SpecializedMove*>(this)->init(ingredients);
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index. Done by the SpecializedMove.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter
 *
 * @tparam <SpecializedMove> name of the specialized move.
 **/
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType>
void MoveForceEquilibriumBase<SpecializedMove>::init(const IngredientsType& ingredients, uint32_t index)
{
  static_cast<SpecializedMove*>(this)->init(ingredients, index);
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index and bond partner. Done by the SpecializedMove.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter
 *
 * @tparam <SpecializedMove> name of the specialized move.
 **/
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType>
void MoveForceEquilibriumBase<SpecializedMove>::init(const IngredientsType& ingredients, uint32_t index, VectorDouble3 dir)
{
  static_cast<SpecializedMove*>(this)->init(ingredients, index, dir);
}

/*****************************************************************************/
/**
 * @brief Check if the move is accepted by the system. Done by the SpecializedMove.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter.
 *
 * @tparam <SpecializedMove> name of the specialized move.
 **/
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType>
void MoveForceEquilibriumBase<SpecializedMove>::check(const IngredientsType& ingredients)
{
  static_cast<SpecializedMove*>(this)->check(ingredients);
}

/*****************************************************************************/
/**
 * @brief Apply the move to the system.
 *
 * @details Here, this is only redirected to the implementation
 * given in the template parameter
 *
 * @tparam <SpecializedMove> name of the specialized move.
 * */
/*****************************************************************************/
template <class SpecializedMove>
template <class IngredientsType>
void MoveForceEquilibriumBase<SpecializedMove>::apply(IngredientsType& ingredients)
{
  static_cast<SpecializedMove*>(this)->apply(ingredients);
}

#endif /*LEMONADE_PM_UPDATER_MOVES_MOVEFORCEEQUILIBRIUMBASE_H*/
