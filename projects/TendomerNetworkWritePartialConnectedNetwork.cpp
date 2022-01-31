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

/****************************************************************************** 
 * based on LeMonADE: https://github.com/LeMonADE-project/LeMonADE/
 * author: Toni Müller
 * email: mueller-toni@ipfdd.de
 * project: LeMonADE-Phantom Modulus
 *****************************************************************************/
#include <iostream>
#include <vector>
#include <bitset>

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFileSubGroup.h>
#include <LeMonADE/feature/FeatureMoleculesIOUnsaveCheck.h>
#include <LeMonADE/feature/FeatureReactiveBonds.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureLabel.h>
#include <LeMonADE/utility/DepthIteratorPredicates.h>
#include <LeMonADE/utility/DepthIterator.h>
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/TaskManager.h>

#include <extern/catchorg/clara/clara.hpp>


#include <LeMonADE_PM/updater/UpdaterReadCrosslinkConnectionsTendomer.h>
#include <LeMonADE_PM/updater/moves/MoveForceEquilibrium.h>
#include <LeMonADE_PM/updater/moves/MoveNonLinearForceEquilibrium.h>
#include <LeMonADE_PM/feature/FeatureCrosslinkConnectionsLookUpTendomers.h>
#include <LeMonADE_PM/analyzer/AnalyzerEquilbratedPosition.h>


int main(int argc, char* argv[]){
	try{
		///////////////////////////////////////////////////////////////////////////////
		///parse options///
		std::string inputBFM("init.bfm");
		std::string outputBFM("Config.bfm.dat");
		std::string inputConnection("BondCreationBreaking.dat");
		bool showHelp = false;
		auto parser
			= clara::detail::Opt(            inputBFM, "inputBFM (=inconfig.bfm)"                        ) ["-i"]["--input"           ] ("(required)Input filename of the bfm file"                                    ).required()
			| clara::detail::Opt(     inputConnection, "inputConnection (=BondCreationBreaking.dat)"     ) ["-d"]["--inputConnection" ] ("used for the time development of the topology. "                             ).required()
			| clara::detail::Opt(           outputBFM, "outputBFM (=Config.bfm)"                         ) ["-o"]["--outputBFM"       ] ("Output filename for the active material of the tendomer network.")
			| clara::Help( showHelp );
		
	    auto result = parser.parse( clara::Args( argc, argv ) );
	    
	    if( !result ) {
	      std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
	      exit(1);
	    }else if(showHelp == true){
	      std::cout << "Analyzer taking a bond table and a table specifying the activity of an object and writes the active part of the tendoemr network."<< std::endl;
          std::cout << "Important: input file for the connections and the active material must fit to each other!!!"<< std::endl;
	      parser.writeToStream(std::cout);
	      exit(0);
	    }else{
	      std::cout << "outputBFM             : " << outputBFM          << std::endl;
	      std::cout << "inputBFM              : " << inputBFM               << std::endl; 
	      std::cout << "inputConnection       : " << inputConnection        << std::endl; 
	    }
		RandomNumberGenerators rng;
		rng.seedAll();
		///////////////////////////////////////////////////////////////////////////////
		///end options parsing
		///////////////////////////////////////////////////////////////////////////////
		//Read in th last Config 
		typedef LOKI_TYPELIST_4(FeatureMoleculesIOUnsaveCheck, FeatureLabel, FeatureReactiveBonds, FeatureAttributes<>) Features;
		typedef ConfigureSystem<VectorInt3,Features, 7> Config;
		typedef Ingredients<Config> Ing;
		Ing myIngredients;
		
		TaskManager taskmanager;
		
		taskmanager.addUpdater( new UpdaterReadBfmFile<Ing>(inputBFM,myIngredients, UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);
		taskmanager.addUpdater( new UpdaterReadCrosslinkConnectionsTendomer<Ing>(myIngredients, inputConnection, 0., 1.0) );
        taskmanager.addAnalyzer( new AnalyzerWriteBfmFile<Ing>(outputBFM, myIngredients, 1) ); 
		//initialize and run
		taskmanager.initialize();
		taskmanager.run(1);
		taskmanager.cleanup();
		std::cout << "Read in conformation and go on to bring it into equilibrium forces..." <<std::endl;
	}
	catch(std::exception& e){
		std::cerr<<"Error:\n"
		<<e.what()<<std::endl;
	}
	catch(...){
		std::cerr<<"Error: unknown exception\n";
	}
	
	return 0;
}
