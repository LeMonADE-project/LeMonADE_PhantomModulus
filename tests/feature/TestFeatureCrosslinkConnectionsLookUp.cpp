/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2021 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
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


/*********************************************************************
 * written by      : Toni Müller
 * email           : mueller-toni@ipfdd.de
 * subprojecttitle : slide ring gels
 *********************************************************************/
#include <iostream>
#include <exception>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/feature/FeatureSystemInformationLinearMeltWithCrosslinker.h>
#include <extern/catch.hpp>

#include <LeMonADE_PM/feature/FeatureCrosslinkConnectionsLookUp.h>

TEST_CASE( "Test class FeatureCrosslinkConnectionsLookUp" ) 
{
    typedef LOKI_TYPELIST_3(FeatureBox, FeatureCrosslinkConnectionsLookUp,FeatureSystemInformationLinearMeltWithCrosslinker ) Features;
    typedef ConfigureSystem<VectorDouble3,Features,4> Config;
    typedef Ingredients<Config> IngredientsType;

  
    std::streambuf* originalBuffer;
    std::ostringstream tempStream;
    //redirect stdout 
    originalBuffer=std::cout.rdbuf();
    std::cout.rdbuf(tempStream.rdbuf());
  
    SECTION(" Test if feature stes up lookup correctly","[FeatureCrosslinkConnectionsLookUp]")
    {
        //setup system 
        IngredientsType ingredients;
        //prepare ingredients
        ingredients.setBoxX(16);
        ingredients.setBoxY(16);
        ingredients.setBoxZ(16);
        ingredients.setPeriodicX(1);
        ingredients.setPeriodicY(1);
        ingredients.setPeriodicZ(1);
        // ingredients.modifyBondset().addBFMclassicBondset();
        //define 
        ingredients.modifyMolecules().addMonomer(6.,6.,6.);
        ingredients.modifyMolecules().addMonomer(6.,4.,6.);
        ingredients.modifyMolecules().addMonomer(6.,8.,6.);
        ingredients.modifyMolecules().addMonomer(4.,6.,6.);
        ingredients.modifyMolecules().addMonomer(8.,6.,6.);
        

        ingredients.modifyMolecules().connect(0,1);
        ingredients.modifyMolecules().connect(0,2);
        ingredients.modifyMolecules().connect(0,3);
        ingredients.modifyMolecules().connect(0,4);

        ingredients.modifyMolecules().addMonomer(6.,4.,6.);
        ingredients.modifyMolecules().addMonomer(6.,4.,6.);
        ingredients.modifyMolecules().connect(1,5);
        ingredients.modifyMolecules().connect(1,6);
        ingredients.modifyMolecules().addMonomer(6.,8.,6.);
        ingredients.modifyMolecules().addMonomer(6.,8.,6.);
        ingredients.modifyMolecules().connect(2,7);
        ingredients.modifyMolecules().connect(2,8);


        ingredients.modifyMolecules().addMonomer(4.,6.,6.);
        ingredients.modifyMolecules().addMonomer(4.,6.,6.);
        ingredients.modifyMolecules().connect(3,9);
        ingredients.modifyMolecules().connect(3,10);
        ingredients.modifyMolecules().addMonomer(8.,6.,6.);
        ingredients.modifyMolecules().addMonomer(8.,6.,6.);
        ingredients.modifyMolecules().connect(4,11);
        ingredients.modifyMolecules().connect(4,12);
        ingredients.setNumOfChains(0);
        ingredients.setNumOfMonomersPerChain(0);
        ingredients.modifyMolecules()[0].setReactive(true); 
        ingredients.modifyMolecules()[0].setNumMaxLinks(4); 
        ingredients.modifyMolecules()[1].setReactive(true); 
        ingredients.modifyMolecules()[1].setNumMaxLinks(4); 
        ingredients.modifyMolecules()[2].setReactive(true); 
        ingredients.modifyMolecules()[2].setNumMaxLinks(4); 
        ingredients.modifyMolecules()[3].setReactive(true); 
        ingredients.modifyMolecules()[3].setNumMaxLinks(4); 
        ingredients.modifyMolecules()[4].setReactive(true); 
        ingredients.modifyMolecules()[4].setNumMaxLinks(4); 



        REQUIRE(ingredients.getMolecules().size()==13 );
        REQUIRE_NOTHROW(ingredients.synchronize(ingredients));
        auto neighbors(ingredients.getCrossLinkNeighborIDs(0));
        REQUIRE(neighbors.size() == 4 ) ;
        REQUIRE(neighbors[0].ID == 1 ) ; 
        REQUIRE(neighbors[1].ID == 2 ) ;
        REQUIRE(neighbors[2].ID == 3 ) ;
        REQUIRE(neighbors[3].ID == 4 ) ;

        REQUIRE(neighbors[0].segDistance == 1 ) ; 
        REQUIRE(neighbors[1].segDistance == 1 ) ; 
        REQUIRE(neighbors[2].segDistance == 1 ) ; 
        REQUIRE(neighbors[3].segDistance == 1 ) ; 
        auto jump(neighbors[0].jump);
        REQUIRE(jump.getX() == Approx(0.) ); 
        REQUIRE(jump.getY() == Approx(0.) ); 
        REQUIRE(jump.getZ() == Approx(0.) ); 

        neighbors=ingredients.getCrossLinkNeighborIDs(1);
        REQUIRE(neighbors[0].ID ==  0 ) ;
        REQUIRE(neighbors.size() == 1 ) ;

        neighbors=ingredients.getCrossLinkNeighborIDs(2);
        REQUIRE(neighbors[0].ID == 0 ) ;
        REQUIRE(neighbors.size() == 1 ) ;

        neighbors=ingredients.getCrossLinkNeighborIDs(3);
        REQUIRE(neighbors[0].ID == 0 ) ;
        REQUIRE(neighbors.size() == 1 ) ;

        neighbors=ingredients.getCrossLinkNeighborIDs(4);
        REQUIRE(neighbors[0].ID == 0 ) ;
        REQUIRE(neighbors.size() == 1 ) ;

        //check if the jumps are correctly calculated
        ingredients.modifyMolecules()[1].modifyVector3D()=( 
            ingredients.getMolecules()[1].getVector3D()+VectorDouble3(16.,64.,-512) );
        REQUIRE_NOTHROW(ingredients.synchronize(ingredients));
        neighbors=ingredients.getCrossLinkNeighborIDs(0);
        jump=(neighbors[0].jump);
        REQUIRE(jump.getX() == Approx(16.) ); 
        REQUIRE(jump.getY() == Approx(64.) ); 
        REQUIRE(jump.getZ() == Approx(-512.) ); 
    }
    //restore cout 
    std::cout.rdbuf(originalBuffer);

}

