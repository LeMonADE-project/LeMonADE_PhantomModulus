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


#ifndef LENONADE_PM_FEATURE_FEATURECROSSLINKCONNECTIONLOOKUPTENDOMERS_H
#define LENONADE_PM_FEATURE_FEATURECROSSLINKCONNECTIONLOOKUPTENDOMERS_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureReactiveBonds.h>
#include <LeMonADE/updater/moves/MoveBase.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE_PM/updater/moves/MoveForceEquilibrium.h>
#include <LeMonADE_PM/utility/neighborX.h>


/*****************************************************************************/
/**
 * @file
 * @date   2021/04/01
 * @author Toni
 *
 * @class FeatureCrosslinkConnectionsLookUpTendomers
 * @brief Creates a lookup for crosslink IDs and their neighbors
 * @details For calculating the equilibrium position of crosslinks in a 
 * a phantom network, one needs to know the crosslink neighbors and the number 
 * of segments between them.  
 * 
 **/
/*****************************************************************************/
///////////////////////////////////////////////////////////////////////////////
//DEFINITION OF THE CLASS TEMPLATE   	                                ///////
//Implementation of the members below					///////
///////////////////////////////////////////////////////////////////////////////

class FeatureCrosslinkConnectionsLookUpTendomers : public Feature {
  
public:
  	//! This Feature requires a monomer_extensions.
	typedef LOKI_TYPELIST_1(MonomerReactivity) monomer_extensions;

	//! check bas connect move - always true 
	template<class IngredientsType>
	bool checkMove(const IngredientsType& ingredients, const MoveBase& move) const { return true;};

	//! synchronize lookup table
	template<class IngredientsType>
	void synchronize(IngredientsType& ingredients) {
	  	fillTables(ingredients);
	};

	//!getter function for the neighboring crosslinks
	std::vector<neighborX> getCrossLinkNeighborIDs(uint32_t CrossLinkID) const{
		if (CrossLinkNeighbors.find(CrossLinkID) == CrossLinkNeighbors.end()){
			std::stringstream errormessage;
			errormessage << "FeatureCrosslinkConnectionsLookUpTendomers::getCrossLinkNeighborIDs Cross Link ID " << CrossLinkID <<" does not exist.";
			throw std::runtime_error(errormessage.str());
		}
		return CrossLinkNeighbors.at(CrossLinkID);
	};

	//! set the jump vector 
	void setCrossLinkNeighborJump(uint32_t CrossLinkID, uint32_t idx, VectorDouble3 vec) {
		if ( idx > CrossLinkNeighbors.at(CrossLinkID).size() ){
			std::stringstream errormessage;
			errormessage << "FeatureCrosslinkConnectionsLookUpTendomers::setCrossLinkNeighborJump neighbor idx  " << idx <<" is to high.";
			throw std::runtime_error(errormessage.str());
		}
		CrossLinkNeighbors.at(CrossLinkID)[idx].jump=vec;
	}

	//!get the ID of crosslinks (determined by nConnections>3 and connected to another crosslink)
	const std::vector<uint32_t>& getCrosslinkIDs() const {return crosslinkIDs;}

private:
  //! convinience function to fill all tables 
  template<class IngredientsType>
  void fillTables(IngredientsType& ingredients);
  //!key - CrossLink ID; value - neighboring cross link and the number of segments to them
  std::map<uint32_t,std::vector< neighborX > > CrossLinkNeighbors;
  //!ID for crosslinks
  std::vector<uint32_t> crosslinkIDs;
};
/**
 *@details  Create look up table 
 **/
template<class IngredientsType>
void FeatureCrosslinkConnectionsLookUpTendomers::fillTables(IngredientsType& ingredients){

	const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
	std::cout << "FeatureCrosslinkConnectionsLookUpTendomers::fillTables" <<std::endl;
	crosslinkIDs.resize(0);
	CrossLinkNeighbors.clear();
    auto nMonomersPerChain(ingredients.getNumMonomersPerChain());
	for (uint32_t i = ingredients.getNumTendomers()*2*nMonomersPerChain;i < molecules.size();i++){
		if( molecules[i].isReactive() && molecules[i].getNumMaxLinks() > 2 ){
			//temporary storage for the neighboring crosslinks
			std::vector<neighborX> NeighborIDs;
			auto posX(molecules[i].getVector3D());
			for (size_t j = 0 ; j < molecules.getNumLinks(i); j++){
				uint32_t tail(i);
				uint32_t head(molecules.getNeighborIdx(i,j));
				VectorDouble3 posHead(molecules[head].getVector3D());
				VectorDouble3 bond(LemonadeDistCalcs::MinImageVector( posX,posHead,ingredients));
				if(bond.getLength() > std::sqrt(10)){
					std::stringstream errormessage;
					errormessage << "FeatureCrosslinkConnectionsLookUpTendomers: Wrong bond " << bond << " between " << i << " and " << head << "\n";
					throw std::runtime_error(errormessage.str());
				}
				VectorDouble3 jumpVector(posHead-bond-posX); // tracks if one bond jumps across periodic images 
				//go to tendomer end
                tail=(head/nMonomersPerChain % 2 ==  0 ) ? head +nMonomersPerChain : head-nMonomersPerChain;
                //cross links are connected by a chain 
                if( molecules.getNumLinks(tail) == 2  ){
                    //get crosslink
                    for (size_t k = 0 ; k < molecules.getNumLinks(tail); k++){
                        uint32_t NextMonomer( molecules.getNeighborIdx(tail,k));
                        if( molecules[NextMonomer].isReactive() && molecules[NextMonomer].getNumMaxLinks() > 2 ){
                            head=NextMonomer; 
                            break;
                        }
                    }
                    posHead=molecules[head].getVector3D();
                    auto posTail=molecules[tail].getVector3D();
                    bond=LemonadeDistCalcs::MinImageVector( posTail,posHead,ingredients);
                    jumpVector+=(posHead-bond-posTail); // tracks if one bond jumps across periodic images 
                    // auto vecJ(jumpVector);
                    // vecJ.setAllCoordinates(
                    //     static_cast<int>(jumpVector.getX()) % 256 , 
                    //     static_cast<int>(jumpVector.getY()) % 256 ,
                    //     static_cast<int>(jumpVector.getZ()) % 256 );

                    // std::cout << "Jumpvector=" << jumpVector << " for ID=" << i << " to " << head  << ": " <<vecJ<<std::endl;
                    auto nSegments(1);
                    NeighborIDs.push_back(neighborX(head, nSegments, jumpVector));

                }
			}
			CrossLinkNeighbors[i]=NeighborIDs;
			crosslinkIDs.push_back(i);
		}
	}
	std::cout << "FeatureCrosslinkConnectionsLookUpTendomers::fillTables.done" <<std::endl; 
}
#endif /*LENONADE_PM_FEATURE_FEATURECROSSLINKCONNECTIONLOOKUP_H*/