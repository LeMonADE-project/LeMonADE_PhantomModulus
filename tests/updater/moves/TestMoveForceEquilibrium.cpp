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
#include <LeMonADE/feature/FeatureSystemInformationLinearMeltWithCrosslinker.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/utility/Vector3D.h>

#include <extern/catch.hpp>

#include <LeMonADE_PM/feature/FeatureCrosslinkConnectionsLookUp.h>
#include <LeMonADE_PM/updater/moves/MoveForceEquilibrium.h>
#include <LeMonADE_PM/feature/FeatureFixedMonomers.h>


TEST_CASE( "Test class MoveForceEquilibrium" ) 
{
    typedef LOKI_TYPELIST_4(FeatureBox, FeatureCrosslinkConnectionsLookUp,FeatureFixedMonomers,FeatureSystemInformationLinearMeltWithCrosslinker) Features;
    typedef ConfigureSystem<VectorDouble3,Features,4> Config;
    typedef Ingredients<Config> IngredientsType;

  
    std::streambuf* originalBuffer;
    std::ostringstream tempStream;
    //redirect stdout 
    originalBuffer=std::cout.rdbuf();
    std::cout.rdbuf(tempStream.rdbuf());
  
    SECTION(" Test if the labels are moved ","[MoveForceEquilibrium]")
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

        ingredients.modifyMolecules()[0].setReactive(true); 
        ingredients.modifyMolecules()[0].setNumMaxLinks(4); 
        for(auto i=1; i < ingredients.getMolecules().size(); i ++){
            ingredients.modifyMolecules()[i].setReactive(true); 
            ingredients.modifyMolecules()[i].setNumMaxLinks(3); 
        }
        for(auto i=2; i < ingredients.getMolecules().size(); i ++)
            ingredients.modifyMolecules()[i].setMovableTag(false);
 


        ingredients.modifyMolecules().addMonomer(4.,6.,6.);
        ingredients.modifyMolecules().addMonomer(4.,6.,6.);
        ingredients.modifyMolecules().connect(3,9);
        ingredients.modifyMolecules().connect(3,10);
        ingredients.modifyMolecules().addMonomer(8.,6.,6.);
        ingredients.modifyMolecules().addMonomer(8.,6.,6.);
        ingredients.modifyMolecules().connect(4,11);
        ingredients.modifyMolecules().connect(4,12);

        REQUIRE(ingredients.getMolecules().size()==13 );
        REQUIRE_NOTHROW(ingredients.synchronize(ingredients));
        //check some basics
        MoveForceEquilibrium move;

        move.init(ingredients,0);
        REQUIRE(move.getShiftVector()==VectorDouble3(0.,0.,0.));
        REQUIRE_NOTHROW(ingredients.synchronize(ingredients));
        ingredients.modifyMolecules()[0].setAllCoordinates(6.,6.,8.);
        move.init(ingredients,0);
        auto vec=move.getShiftVector();
        REQUIRE(vec.getY() == 0 );
        REQUIRE(vec.getX() == 0 );
        REQUIRE(vec.getZ() == Approx(-2.) );
        ingredients.modifyMolecules()[0].setAllCoordinates(12.0,17.,3.0);
        move.init(ingredients,0);
        move.check(ingredients);
        move.apply(ingredients);
        vec=move.getShiftVector();
        REQUIRE(vec.getX() == Approx( -6.));
        REQUIRE(vec.getY() == Approx(-11.));
        REQUIRE(vec.getZ() == Approx( 3.));
    }
    //restore cout 
    std::cout.rdbuf(originalBuffer);

}

