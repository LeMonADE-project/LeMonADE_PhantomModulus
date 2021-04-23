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
#include <LeMonADE/utility/Vector3D.h>

#include <extern/catch.hpp>


#include <LeMonADE_PM/updater/moves/MoveNonLinearForceEquilibrium.h>


TEST_CASE( "Test class MoveNonLinearForceEquilibrium" ) 
{
 
    std::streambuf* originalBuffer;
    std::ostringstream tempStream;
    //redirect stdout 
    originalBuffer=std::cout.rdbuf();
    std::cout.rdbuf(tempStream.rdbuf());
  
    SECTION("Check setter/getter in  data file ","[MoveNonLinearForceEquilibrium_GETSET]")
    {

        MoveNonLinearForceEquilibrium move; 
        move.setFilename("MaxisCruve.dat");
        REQUIRE(move.getFilename() == "MaxisCruve.dat");
        move.setRelaxationParameter(123.28);
        REQUIRE(move.getRelaxationParameter() == 123.28);
        move.setRelaxationParameter(16.);
        REQUIRE(19.1530666667 == Approx(move.FE(VectorDouble3(0.5, 0.,0.)).getLength()) );
    }
    SECTION ("Check the reading of a file1", "[MoveNonLinearForceEquilibrium_READIN1]")
    {
        //read in a perfectly spaces curve
        std::string filename("TomsCurve.dat");
        std::ofstream out(filename);
        const double N=32.;
        const double b=2.68;
        for(auto i=0; i < 100; i++ ){
            double force(static_cast<double>(i)*3./(N*b*b));
            //force extension :
            out << force  << "\t"<<i << "\n"; 
        }
        out.close();
        MoveNonLinearForceEquilibrium move(filename); 
        REQUIRE(move.EF(VectorDouble3( 0.0,0.,0.)).getLength()==Approx(0.));
        REQUIRE(move.EF(VectorDouble3( 5.0,0.,0.)).getLength()==Approx(0.06526370015));
        REQUIRE(move.EF(VectorDouble3( 5.4,0.,0.)).getLength()==Approx(0.06526370015));
        REQUIRE(move.EF(VectorDouble3( 4.6,0.,0.)).getLength()==Approx(0.06526370015));
        REQUIRE(0==remove(filename.c_str()));    
    }
    SECTION ("Check the reading of a file2", "[MoveNonLinearForceEquilibrium_READIN2]")
    {
        //read in a perfectly spaces curve
        std::string filename("SmotsCurve.dat");
        std::ofstream out(filename);
        const double N=32.;
        const double b=2.68;
        for(auto i=99; i >=0; i-- ){
            double force(static_cast<double>(i)*3./(N*b*b));
            //force extension :
            out << force  << "\t"<<i << "\n"; 
        }
        out.close();
        MoveNonLinearForceEquilibrium move(filename); 
        REQUIRE(move.EF(VectorDouble3( 0.0,0.,0.)).getLength()==Approx(0.));
        REQUIRE(move.EF(VectorDouble3( 5.0,0.,0.)).getLength()==Approx(0.06526370015));
        REQUIRE(move.EF(VectorDouble3( 5.4,0.,0.)).getLength()==Approx(0.06526370015));
        REQUIRE(move.EF(VectorDouble3( 4.6,0.,0.)).getLength()==Approx(0.06526370015));
        REQUIRE(0==remove(filename.c_str()));    
    }
    SECTION ("Check the reading of a file3", "[MoveNonLinearForceEquilibrium_READIN3]")
    {
        //read in  curve with bigger spaces 
        std::string filename("TmCurve.dat");
        std::ofstream out(filename);
        const double N=32.;
        const double b=2.68;
        for(auto i=0; i <100; i+=5 ){
            double force(static_cast<double>(i)*3./(N*b*b));
            //force extension :
            out << force  << "\t"<<i << "\n"; 
        }
        out.close();
        MoveNonLinearForceEquilibrium move(filename); 
        REQUIRE(move.EF(VectorDouble3( 0.0,0.,0.)).getLength()==Approx(0.));
        REQUIRE(move.EF(VectorDouble3( 5.0,0.,0.)).getLength()==Approx(0.06526370015));
        REQUIRE(move.EF(VectorDouble3( 5.4,0.,0.)).getLength()==Approx(0.06526370015));
        REQUIRE(move.EF(VectorDouble3( 4.6,0.,0.)).getLength()==Approx(0.06526370015));
        REQUIRE(0==remove(filename.c_str()));    
    }

    //restore cout 
    std::cout.rdbuf(originalBuffer);

}

