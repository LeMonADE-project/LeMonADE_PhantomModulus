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
 * subprojecttitle : Phantom modulus
 *********************************************************************/
#include <iostream>
#include <exception>

#include <LeMonADE/core/Molecules.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureSystemInformationLinearMeltWithCrosslinker.h>

#include <LeMonADE/utility/Vector3D.h>

#include <extern/catch.hpp>

#include <LeMonADE_PM/updater/UpdaterReadCrosslinkConnections.h>
#include <LeMonADE_PM/feature/FeatureCrosslinkConnectionsLookUp.h>



TEST_CASE( "Test class UpdaterReadCrosslinkConnections" ) 
{
    typedef LOKI_TYPELIST_3(FeatureBox, FeatureSystemInformationLinearMeltWithCrosslinker,FeatureCrosslinkConnectionsLookUp) Features;
    typedef ConfigureSystem<VectorDouble3,Features,4> Config;
    typedef Ingredients<Config> IngredientsType;

  
    std::streambuf* originalBuffer;
    std::ostringstream tempStream;
    //redirect stdout 
    originalBuffer=std::cout.rdbuf();
    std::cout.rdbuf(tempStream.rdbuf());
  
    SECTION(" Test if the labels are moved ","[UpdaterReadCrosslinkConnections]")
    {
        //prepare input file 
        const std::string filename("bondTable.dat");
        std::ofstream out(filename); 
        out <<"# \n";
        out <<"\n";
        //   Time >>  ChainID >>    MonID1 >>       P1X >>     P1Y >>     P1Z >>   MonID2 >>      P2X >>     P2Y >>     P2Z
        out << 17 << " " << 1 << " " << 12 << " " << 12 << " "<< 3 << " "<< 2 << " "<<  0 << " "<< 12 << " "<< 1 << " "<< 3 <<"\n";
        out << 17 << " " << 1 << " " << 13 << " " << 12 << " "<< 3 << " "<< 2 << " "<< 12 << " "<< 12 << " "<< 1 << " "<< 3 <<"\n";

        out << 19 << " " << 2 << " " << 14 << " " << 12 << " "<< 3 << " "<< 2 << " "<< 0 << " "<< 0 << " "<< 4 << " "<< 3 <<"\n";
        out << 19 << " " << 2 << " " << 12 << " " << 12 << " "<< 3 << " "<< 2 << " "<< 14 << " "<< 12 << " "<< 4 << " "<< 3 <<"\n";

        out << 21 << " " << 3 << " " << 12 << " " <<  2 << " "<< 3 << " "<< 2 << " "<< 0 << " "<< 12 << " "<< 4 << " "<< 7 <<"\n";
        out << 23 << " " << 3 << " " << 15 << " " <<  2 << " "<< 3 << " "<< 2 << " "<< 12 << " "<< 12 << " "<< 4 << " "<< 7 <<"\n";

        out << 27 << " " << 4 << " " << 16 << " " << 45 << " "<< 3 << " "<< 2 << " "<< 0 << " "<<  5 << " "<< 4 << " "<< 9 <<"\n";
        out << 29 << " " << 4 << " " << 12 << " " << 45 << " "<< 3 << " "<< 2 << " "<< 16 << " "<<  5 << " "<< 4 << " "<< 9 <<"\n";

        out.close();
        //setup system 
        IngredientsType ingredients;
        //prepare ingredients
        ingredients.setBoxX(16);
        ingredients.setBoxY(16);
        ingredients.setBoxZ(16);
        ingredients.setPeriodicX(1);
        ingredients.setPeriodicY(1);
        ingredients.setPeriodicZ(1);
        ingredients.setNumOfChains(12);
        ingredients.setNumOfCrosslinks(5);
        ingredients.setFunctionality(4);
        ingredients.setNumOfMonomersPerChain(1);
        ingredients.setNumOfMonomersPerCrosslink(1);
        //define 
        //chains 
        ingredients.modifyMolecules().addMonomer(6.,5.,6.);//0
        ingredients.modifyMolecules().addMonomer(6.,7.,6.);//1
        ingredients.modifyMolecules().addMonomer(5.,6.,6.);//2
        ingredients.modifyMolecules().addMonomer(7.,6.,6.);//3

        ingredients.modifyMolecules().addMonomer(6.,4.,6.);//4
        ingredients.modifyMolecules().addMonomer(6.,4.,6.);//5
        ingredients.modifyMolecules().addMonomer(6.,8.,6.);//6
        ingredients.modifyMolecules().addMonomer(6.,8.,6.);//7
        ingredients.modifyMolecules().addMonomer(4.,6.,6.);//8
        ingredients.modifyMolecules().addMonomer(4.,6.,6.);//9
        ingredients.modifyMolecules().addMonomer(8.,6.,6.);//10
        ingredients.modifyMolecules().addMonomer(8.,6.,6.);//11

        //crosslinks
        ingredients.modifyMolecules().addMonomer(6.,6.,6.);//12
        ingredients.modifyMolecules().addMonomer(6.,4.,6.);//13
        ingredients.modifyMolecules().addMonomer(6.,8.,6.);//14
        ingredients.modifyMolecules().addMonomer(4.,6.,6.);//15
        ingredients.modifyMolecules().addMonomer(8.,6.,6.);//16
        
        ingredients.modifyMolecules().connect(12,0);
        ingredients.modifyMolecules().connect(12,1);
        ingredients.modifyMolecules().connect(12,2);
        ingredients.modifyMolecules().connect(12,3);
        ingredients.modifyMolecules().connect(13,0);
        ingredients.modifyMolecules().connect(14,1);
        ingredients.modifyMolecules().connect(15,2);
        ingredients.modifyMolecules().connect(16,3);


        ingredients.modifyMolecules().connect(13,4);
        ingredients.modifyMolecules().connect(13,5);
        ingredients.modifyMolecules().connect(14,6);
        ingredients.modifyMolecules().connect(14,7);
        ingredients.modifyMolecules().connect(15,8);
        ingredients.modifyMolecules().connect(15,9);
        ingredients.modifyMolecules().connect(16,10);
        ingredients.modifyMolecules().connect(16,11);
        
        // for (auto i=0; i < ingredients.getMolecules().size(); i++){
        for (auto i=0; i < 4; i++){
            ingredients.modifyMolecules()[i].setReactive(true); 
            ingredients.modifyMolecules()[i].setNumMaxLinks(2); 
        }

        ingredients.modifyMolecules()[12].setReactive(true); 
        ingredients.modifyMolecules()[12].setNumMaxLinks(4); 
        ingredients.modifyMolecules()[13].setReactive(true); 
        ingredients.modifyMolecules()[13].setNumMaxLinks(4); 
        ingredients.modifyMolecules()[14].setReactive(true); 
        ingredients.modifyMolecules()[14].setNumMaxLinks(4); 
        ingredients.modifyMolecules()[15].setReactive(true); 
        ingredients.modifyMolecules()[15].setNumMaxLinks(4); 
        ingredients.modifyMolecules()[16].setReactive(true); 
        ingredients.modifyMolecules()[16].setNumMaxLinks(4); 

        REQUIRE(ingredients.getMolecules().size()==17 );
        REQUIRE_NOTHROW(ingredients.synchronize(ingredients));
        UpdaterReadCrosslinkConnections<IngredientsType> updater(ingredients, filename, 1., 0.00);
        updater.initialize();
        for(auto i =0; i < 10 ; i++)
            updater.execute();
        REQUIRE(ingredients.getMolecules().areConnected(12,0) );
        for(auto i =0; i < 10 ; i++)
            updater.execute();
        REQUIRE(ingredients.getMolecules().areConnected(0,13) );
        
        //read whole file 
        while(updater.execute()){}
        //check connections
        REQUIRE(ingredients.getMolecules().areConnected(12,1) );
        REQUIRE(ingredients.getMolecules().areConnected(12,2) );
        REQUIRE(ingredients.getMolecules().areConnected(12,3) );

        REQUIRE(ingredients.getMolecules().areConnected(14,1) );
        REQUIRE(ingredients.getMolecules().areConnected(15,2) );
        REQUIRE(ingredients.getMolecules().areConnected(16,3) );
        
        REQUIRE(0==remove(filename.c_str()));    
    }
    //restore cout 
    std::cout.rdbuf(originalBuffer);

}

