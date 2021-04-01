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


#ifndef LENONADE_PM_FEATURE_FEATURECROSSLINKCONNECTIONLOOKUP_H
#define LENONADE_PM_FEATURE_FEATURECROSSLINKCONNECTIONLOOKUP_H

#include <LeMonADE/feature/Feature.h>
#include <LeMonADE/feature/FeatureConnectionSc.h>
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
 * @class FeatureCrosslinkConnectionsLookUp
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

class FeatureCrosslinkConnectionsLookUp : public Feature {
  
public:

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
		if (CrossLinkNeighbors.find(CrossLinkID) == CrossLinkNeighbors.end())
		{
			std::stringstream errormessage;
			errormessage << "FeatureCrosslinkConnectionsLookUp::getCrossLinkNeighborIDs Cross Link ID " << CrossLinkID <<" does not exist.";
			throw std::runtime_error(errormessage.str());
		}
		return CrossLinkNeighbors.at(CrossLinkID);
	};

	// std::vector<neighborX >& modifyCrossLinkNeighborIDs(uint32_t CrossLinkID) const{
	// 	return CrossLinkNeighbors[CrossLinkID];
	// };

private:
  //! convinience function to fill all tables 
  template<class IngredientsType>
  void fillTables(IngredientsType& ingredients);
  //!key - CrossLink ID; value - neighboring cross link and the number of segments to them
  std::map<uint32_t,std::vector< neighborX > > CrossLinkNeighbors;
};
/**
 *@details  Create look up table 
 **/
template<class IngredientsType>
void FeatureCrosslinkConnectionsLookUp::fillTables(IngredientsType& ingredients){
	const typename IngredientsType::molecules_type& molecules=ingredients.getMolecules();
	std::cout << "FeatureCrosslinkConnectionsLookUp::fillTables" <<std::endl;
	for (uint32_t i = 0 ;i < molecules.size();i++){
			//find next crosslink
			if( molecules.getNumLinks(i) > 2 ){
				//temporary storage for the neighboring crosslinks
				std::vector<neighborX> NeighborIDs;
				auto posX(molecules[i].getVector3D());
				for (size_t j = 0 ; j < molecules.getNumLinks(i); j++){
					uint32_t tail(i);
					uint32_t head(molecules.getNeighborIdx(i,j));
					auto posHead(molecules[head].getVector3D());
					auto bond(LemonadeDistCalcs::MinImageVector( posX,posHead,ingredients));
					auto jumpVector(posHead-bond-posX); // tracks if one bond jumps across periodic images 
					//direct connection of two cross links
					if (molecules.getNumLinks(head) > 2) {
						NeighborIDs.push_back( neighborX(head, 1, jumpVector) );
					}else{ 
						uint32_t nSegments(1);
						//cross links are connected by a chain 
						while( molecules.getNumLinks(head) == 2 ){
							//find next head 
							for (size_t k = 0 ; k < molecules.getNumLinks(head); k++){
								uint32_t NextMonomer( molecules.getNeighborIdx(head,k));
								if ( NextMonomer != tail ) {
									tail=head;
									head=NextMonomer; 
									break;
								}
							}
							posHead=molecules[head].getVector3D();
							bond=LemonadeDistCalcs::MinImageVector( posX,posHead,ingredients);
							jumpVector+=(posHead-bond); // tracks if one bond jumps across periodic images 
							nSegments++;
							//a cross link has more than 2 connections
							if (molecules.getNumLinks(head) > 2) {
								NeighborIDs.push_back(neighborX(head, nSegments, jumpVector));
								break;
							}
						}
                    }	
				}
				CrossLinkNeighbors[i]=NeighborIDs;
			}
	}
	std::cout << "FeatureCrosslinkConnectionsLookUp::fillTables.done" <<std::endl; 
}
#endif /*LENONADE_PM_FEATURE_FEATURECROSSLINKCONNECTIONLOOKUP_H*/