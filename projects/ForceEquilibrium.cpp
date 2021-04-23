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
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/feature/FeatureMoleculesIOUnsaveCheck.h>
#include <LeMonADE/feature/FeatureReactiveBonds.h>
#include <LeMonADE/feature/FeatureBox.h>
#include <LeMonADE/feature/FeatureSystemInformationLinearMeltWithCrosslinker.h>
#include <LeMonADE/utility/TaskManager.h>

#include <extern/catchorg/clara/clara.hpp>

#include <LeMonADE_PM/updater/UpdaterForceBalancedPosition.h>
#include <LeMonADE_PM/updater/UpdaterReadCrosslinkConnections.h>
#include <LeMonADE_PM/updater/moves/MoveForceEquilibrium.h>
#include <LeMonADE_PM/updater/moves/MoveNonLinearForceEquilibrium.h>
#include <LeMonADE_PM/feature/FeatureCrosslinkConnectionsLookUp.h>
#include <LeMonADE_PM/analyzer/AnalyzerEquilbratedPosition.h>


int main(int argc, char* argv[]){
	try{
		///////////////////////////////////////////////////////////////////////////////
		///parse options///
		std::string inputBFM("init.bfm");
		std::string outputDataPos("CrosslinkPosition.dat");
		std::string outputDataDist("ChainExtensionDistribution.dat");
		std::string inputConnection("BondCreationBreaking.dat");
		std::string feCurve("");
		double relaxationParameter(10.);
		double threshold(0.5);
		double stepwidth(1.0);
		double minConversion(50.0);
		bool custom(false);
		
		bool showHelp = false;
		auto parser
			= clara::detail::Opt(            inputBFM, "inputBFM (=inconfig.bfm)"                        ) ["-i"]["--input"          ] ("(required)Input filename of the bfm file"                                    ).required()
			| clara::detail::Opt(     inputConnection, "inputConnection (=BondCreationBreaking.dat)"     ) ["-d"]["--inputConnection"] ("used for the time development of the topology. "                             ).required()
			| clara::detail::Opt(       outputDataPos, "outputDataPos (=CrosslinkPosition.dat)"          ) ["-o"]["--outputPos"      ] ("(optional) Output filename of the crosslink ID and the equilibrium Position.").optional()
			| clara::detail::Opt(      outputDataDist, "outputDataDist (=ChainExtensionDistribution.dat)") ["-c"]["--outputDist"     ] ("(optional) Output filename of the chain extension distribution."             ).optional()
			| clara::detail::Opt(           stepwidth, "stepwidth"                                       ) ["-s"]["--stepwidth"      ] ("(optional) Width for the increase in percentage. Default: 1%."               ).optional()
			| clara::detail::Opt(       minConversion, "minConversion"                                   ) ["-u"]["--minConversion"  ] ("(optional) Minimum conversion to be read in. Default: 50%."                  ).optional()
			| clara::detail::Opt(           threshold, "threshold"                                       ) ["-t"]["--threshold"      ] ("(optional) Threshold of the average shift. Default 0.5 ."                    ).optional()
			| clara::detail::Opt(             feCurve, "feCurve (="")"                                   ) ["-f"]["--feCurve"        ] ("(optional) Force-Extension curve. Default \"\"."                             ).optional()
			| clara::detail::Opt( relaxationParameter, "relaxationParameter (=10)"                       ) ["-r"]["--relax"          ] ("(optional) Relaxation parameter. Default 10.0 ."                             ).optional()
			| clara::Help( showHelp );
		
	    auto result = parser.parse( clara::Args( argc, argv ) );
	    
	    if( !result ) {
	      std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
	      exit(1);
	    }else if(showHelp == true){
	      std::cout << "Simulator to connect linear chains with single monomers of certain functionality"<< std::endl;
	      parser.writeToStream(std::cout);
	      exit(0);
	    }else{
	      std::cout << "outputData            : " << outputDataPos          << std::endl;
	      std::cout << "outputDataDist        : " << outputDataDist         << std::endl;
	      std::cout << "inputBFM              : " << inputBFM               << std::endl; 
	      std::cout << "inputConnection       : " << inputConnection        << std::endl; 
	      std::cout << "stepwidth             : " << stepwidth              << std::endl;
	      std::cout << "minConversion         : " << minConversion          << std::endl;
	      std::cout << "threshold             : " << threshold              << std::endl; 
		  std::cout << "feCurve               : " << feCurve                << std::endl;
	    }
		if (! feCurve.empty()) custom=true;
		
		///////////////////////////////////////////////////////////////////////////////
		///end options parsing
		///////////////////////////////////////////////////////////////////////////////
		//Read in th last Config 
		typedef LOKI_TYPELIST_3(FeatureMoleculesIOUnsaveCheck, FeatureSystemInformationLinearMeltWithCrosslinker, FeatureReactiveBonds) Features;
		typedef ConfigureSystem<VectorInt3,Features, 7> Config;
		typedef Ingredients<Config> Ing;
		Ing myIngredients;
		
		TaskManager taskmanager;
		
		taskmanager.addUpdater( new UpdaterReadBfmFile<Ing>(inputBFM,myIngredients, UpdaterReadBfmFile<Ing>::READ_LAST_CONFIG_SAVE),0);

		//initialize and run
		taskmanager.initialize();
		taskmanager.run(1);
		taskmanager.cleanup();
		std::cout << "Read in conformation and go on to bring it into equilibrium forces..." <<std::endl;
		//the foce equilibrium is reached off lattice ( no integer values for the positions )
		typedef LOKI_TYPELIST_3(FeatureBox, FeatureCrosslinkConnectionsLookUp ,FeatureSystemInformationLinearMeltWithCrosslinker) Features2;
		typedef ConfigureSystem<VectorDouble3,Features2, 7> Config2;
		typedef Ingredients<Config2> Ing2;
		Ing2 myIngredients2;
		
		myIngredients2.setBoxX(myIngredients.getBoxX());
		myIngredients2.setBoxY(myIngredients.getBoxY());
		myIngredients2.setBoxZ(myIngredients.getBoxZ());
		myIngredients2.setPeriodicX(myIngredients.isPeriodicX());
		myIngredients2.setPeriodicY(myIngredients.isPeriodicY());
		myIngredients2.setPeriodicZ(myIngredients.isPeriodicZ());
		myIngredients2.modifyMolecules().resize(myIngredients.getMolecules().size());
		myIngredients2.modifyMolecules().setAge(myIngredients.getMolecules().getAge());
		myIngredients2.setNumOfChains              (myIngredients.getNumOfChains());
		myIngredients2.setNumOfCrosslinks          (myIngredients.getNumOfCrosslinks());
		myIngredients2.setNumOfMonomersPerChain    (myIngredients.getNumOfMonomersPerChain());
		myIngredients2.setNumOfMonomersPerCrosslink(myIngredients.getNumOfMonomersPerCrosslink());
		myIngredients2.setFunctionality            (myIngredients.getFunctionality());
		
		for(size_t i = 0; i< myIngredients.getMolecules().size();i++){
			myIngredients2.modifyMolecules()[i].modifyVector3D()=myIngredients.getMolecules()[i].getVector3D();
			myIngredients2.modifyMolecules()[i].setReactive(myIngredients.getMolecules()[i].isReactive());
			myIngredients2.modifyMolecules()[i].setNumMaxLinks(myIngredients.getMolecules()[i].getNumMaxLinks());
			for (size_t j = 0 ; j < myIngredients.getMolecules().getNumLinks(i);j++){
				uint32_t neighbor(myIngredients.getMolecules().getNeighborIdx(i,j));
				if( ! myIngredients2.getMolecules().areConnected(i,neighbor) )
					myIngredients2.modifyMolecules().connect(i,neighbor);
			}
		}
		myIngredients2.synchronize();

		TaskManager taskmanager2;
		//read bonds and positions stepwise
		taskmanager2.addUpdater( new UpdaterReadCrosslinkConnections<Ing2>(myIngredients2, inputConnection, stepwidth, minConversion) );
		if (custom) 
			taskmanager2.addUpdater( new UpdaterForceBalancedPosition<Ing2,MoveForceEquilibrium>(myIngredients2, threshold) );
		else {
			auto updater = new UpdaterForceBalancedPosition<Ing2,MoveNonLinearForceEquilibrium>(myIngredients2, threshold) ;
			updater->setFilename(feCurve);
			updater->setRelaxationParameter(relaxationParameter);
			taskmanager2.addUpdater( updater );
		}
		taskmanager2.addAnalyzer(new AnalyzerEquilbratedPosition<Ing2>(myIngredients2,outputDataPos, outputDataDist));
		//initialize and run
		taskmanager2.initialize();
		taskmanager2.run();
		taskmanager2.cleanup();
		
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
