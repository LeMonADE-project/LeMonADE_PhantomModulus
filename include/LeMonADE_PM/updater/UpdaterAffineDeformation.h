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
#ifndef LEMONADE_PM_UPDATER_UPDATERAFFINEDEFORMATION_H
#define LEMONADE_PM_UPDATER_UPDATERAFFINEDEFORMATION_H


#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE_PM/utility/neighborX.h>
#include <vector>
#include <cmath>
 /**
 * @class UpdaterAffineDeformation
 * @brief This class performs a affine deformation accordign to the stretching factor (lambda)
 * for the positions and the jump vector. 
 * @tparam IngredientsType
 */

 template <class IngredientsType >
class UpdaterAffineDeformation:public AbstractUpdater
{
public:
    //! constructor for UpdaterAffineDeformation
    UpdaterAffineDeformation(IngredientsType& ing_, double stretching_factor_, double prestrain_factorX_=1, double prestrain_factorY_=1, double prestrain_factorZ_=1 ):
    ing(ing_),
    stretching_factor(stretching_factor_),
    stretching_factor_XY(1./std::sqrt(stretching_factor)),
    prestrain_factorX(prestrain_factorX_),
    prestrain_factorY(prestrain_factorY_),
    prestrain_factorZ(prestrain_factorZ_)
    {
        std::cout << "UpdaterAffineDeformation::constructor:\n"
                  << "stretching factor =       " << stretching_factor    << "\n"
                  << "stretching factor in yz = " << stretching_factor_XY << "\n"
                  << "prestrain factor in x = " << prestrain_factorX << "\n"
                  << "prestrain factor in y = " << prestrain_factorY << "\n"
                  << "prestrain factor in z = " << prestrain_factorZ << "\n";
    };
    
    virtual void initialize(){};
    virtual bool execute();
    virtual void cleanup(){};  

    //! setter function for the stretching factor 
    void setStretchingFactor(double stretching_factor_ ) {
        stretching_factor=stretching_factor_;
        stretching_factor_XY = 1./std::sqrt(stretching_factor) ; 
    }

private:
    //!copy of the main container for the system informations 
    IngredientsType& ing;
    
    //! threshold for the certainty 
    double stretching_factor;
    double stretching_factor_XY;
    //! factors for a prestrain or a correction factor for the elastic inbalance in the spatial direction 
    double prestrain_factorX,prestrain_factorY,prestrain_factorZ;
    VectorDouble3 deform( VectorDouble3 vec ){

        return VectorDouble3(vec.getX()*stretching_factor   *prestrain_factorX,
                             vec.getY()*stretching_factor_XY*prestrain_factorY,
                             vec.getZ()*stretching_factor_XY*prestrain_factorZ);
    }
};
template <class IngredientsType>
bool UpdaterAffineDeformation<IngredientsType>::execute(){
    std::cout << "UpdaterAffineDeformation<IngredientsType>::initialize():"<< std::endl;
    //adjusting the box size is not neccessary, because it is used only once in the FeatureCrosslinkConnections*
    //there the jump vectors are calculated
    
    //The jump vector needs to be adjusted according the deformation labmda! 
    auto CrossLinkIDs(ing.getCrosslinkIDs());
    for (uint32_t i =0 ; i<CrossLinkIDs.size() ; i++){
        uint32_t ID(CrossLinkIDs[i]);
        std::vector<neighborX> Neighbors(ing.getCrossLinkNeighborIDs(ID) );
        int32_t number_of_neighbors(Neighbors.size());
        if (number_of_neighbors > 0) {
            for (size_t j = 0; j < number_of_neighbors; j++){
                if (i < 20 ) 
                    std::cout << "ID= "<< i << " initial jump " << Neighbors[j].jump << " ";  
                ing.setCrossLinkNeighborJump(ID,j,deform(Neighbors[j].jump));   
                if (i < 20 ) 
                    std::cout << "ID= "<< i << " deformed jump " << ing.getCrossLinkNeighborIDs(ID)[j].jump << "\n";   
            }
        }
    }

    //The initial posistions of all monomers needs to be adjusted according the deformation labmda! 
    for(size_t i = 0; i< ing.getMolecules().size();i++){
        if (i < 20 ) 
            std::cout << "ID="<< i << "initial position " << ing.getMolecules()[i].getVector3D() << " ";  
        ing.modifyMolecules()[i].modifyVector3D()=deform(ing.getMolecules()[i].getVector3D());
        if (i < 20 ) 
            std::cout << " deformed position " << ing.getMolecules()[i].getVector3D() << "\n";  
    }
    std::cout << "UpdaterAffineDeformation<IngredientsType>::initialize():done.\n";
}
#endif /*LEMONADE_PM_UPDATER_UPDATERAFFINEDEFORMATION_H*/