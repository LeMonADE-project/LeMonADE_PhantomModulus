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
        std::string activeComponent("active_00001.dat");		
		bool showHelp = false;
		auto parser
			= clara::detail::Opt(            inputBFM, "inputBFM (=inconfig.bfm)"                        ) ["-i"]["--input"           ] ("(required)Input filename of the bfm file"                                    ).required()
			| clara::detail::Opt(     inputConnection, "inputConnection (=BondCreationBreaking.dat)"     ) ["-d"]["--inputConnection" ] ("used for the time development of the topology. "                             ).required()
            | clara::detail::Opt(     activeComponent, "activeComponents (=active_00001.dat)"            ) ["-a"]["--activeComponents"] ("sets the active components . "                                              ).required()
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
          std::cout << "activeComponent       : " << activeComponent        << std::endl; 
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

		//initialize and run
		taskmanager.initialize();
		taskmanager.run(1);
		taskmanager.cleanup();
		std::cout << "Read in conformation and go on to bring it into equilibrium forces..." <<std::endl;
        // activeObjects are first all cross links and afterwards all chain: 0-sol, 1-gel/pending, 2-active 
        std::ifstream in(activeComponent);
        std::string line;

		auto nCrossLinks(myIngredients.getNumCrossLinkers());
		auto nChains(myIngredients.getNumTendomers());
		auto nSegments(myIngredients.getNumMonomersPerChain()*2);
		auto nChainMonomers(nChains*nSegments);
		auto objectID(0);
		auto nActiveCrossLinks(0);
		auto nActiveTendomers(0);
		std::cout <<"Filestart:\n" ;
        while (in.good() &&  ! in.eof() ){
            getline(in,line);
            std::stringstream ss;
            uint32_t activeTag;
            ss << line;
			// std::cout << line <<"\n";
			if (! line.empty()) {
				ss >>activeTag;
				auto tmpAttribute(0);
				if(objectID < nCrossLinks ){
					auto IdX(nChainMonomers+objectID);
					#ifdef DEBUG
						std::cout <<  objectID << " "  << activeTag <<  " " << IdX <<  " " << IdX << std::endl;
					#endif
					if (activeTag == 2 || activeTag == 1 ){
						tmpAttribute=1;	
						nActiveCrossLinks++;
					}
					myIngredients.modifyMolecules()[nChainMonomers+objectID].setAttributeTag(tmpAttribute);
				}else{
					auto IdCStart((objectID-nCrossLinks)*nSegments);
					auto IdCEnd( (objectID-nCrossLinks+1)*nSegments-1);
					#ifdef DEBUG
						std::cout << objectID << " "  << activeTag <<" " << IdCStart << " " << IdCEnd << std::endl;
					#endif
					if (activeTag == 2 || activeTag == 1 ) {
						nActiveTendomers++;
						tmpAttribute=1;
					}
					for (auto i =IdCStart ; i <= IdCEnd; i++ )
						myIngredients.modifyMolecules()[i].setAttributeTag(tmpAttribute);
				}
				objectID++;
			}
        }
		std::cout <<"Fileend:\n" ;
		// for()
		myIngredients.setNumTendomers           (nActiveTendomers);
		myIngredients.setNumCrossLinkers        (nActiveCrossLinks);
		myIngredients.setNumMonomersPerChain    (myIngredients.getNumMonomersPerChain());
		myIngredients.setNumLabelsPerTendomerArm(myIngredients.getNumLabelsPerTendomerArm());
		myIngredients.synchronize();
		

		MonomerGroup<typename Ing::molecules_type> subgroup(myIngredients.getMolecules());
		subgroup.clear();
		hasThisType<1> predicate;
		for(size_t n=0;n<myIngredients.getMolecules().size();n++)
		{
			if(predicate(myIngredients.getMolecules(),n)==true) subgroup.push_back(n);
		}
		myIngredients.modifyMolecules()=subgroup.copyGroup();

		for(size_t n=0;n<myIngredients.getMolecules().size();n++)
		{
			myIngredients.modifyMolecules()[n].setReactive(false);
			myIngredients.modifyMolecules()[n].setNumMaxLinks(0);
			if( n < nActiveTendomers*nSegments ){
				if (n%(nSegments/2) == 0 ) {
					myIngredients.modifyMolecules()[n].setReactive(true);
					myIngredients.modifyMolecules()[n].setNumMaxLinks(2);
				}
			}else {
				myIngredients.modifyMolecules()[n].setReactive(true);
				myIngredients.modifyMolecules()[n].setNumMaxLinks(4);
			}
		}		

		myIngredients.synchronize();
		// auto BoxX(myIngredients.getBoxX());
		// auto BoxY(myIngredients.getBoxY());
		// auto BoxZ(myIngredients.getBoxZ());
		// for (auto i=nActiveTendomers*2*myIngredients.getNumMonomersPerChain(); i < myIngredients.getMolecules().size(); i++ ){
		// 	auto vec(myIngredients.getMolecules()[i].getVector3D());
		// 	auto new_vec(vec);
		// 	new_vec.setAllCoordinates( LemonadeDistCalcs::fold(vec.getX(),BoxX), 
		// 		 					   LemonadeDistCalcs::fold(vec.getX(),BoxY),
		// 							   LemonadeDistCalcs::fold(vec.getX(),BoxZ));
		// 	myIngredients.modifyMolecules()[i].modifyVector3D() = new_vec;
		// }
		// myIngredients.synchronize();

		TaskManager taskmanager2;
		taskmanager2.addAnalyzer( new AnalyzerWriteBfmFile<Ing >(outputBFM, myIngredients, 1) ); 
		//initialize and run
		if ( nActiveTendomers > 0  ){
			taskmanager2.initialize();
			taskmanager2.run(1);
			taskmanager2.cleanup();
		}
		// AnalyzerWriteBfmFileSubGroup<Ing,hasThisType<1> >  writer(outputBFM, myIngredients, 1, hasThisType<1>());
		// writer.initialize();
		
		// TaskManager taskmanager2;
		// taskmanager2.addAnalyzer( new AnalyzerWriteBfmFileSubGroup<Ing,hasThisType<1> >(outputBFM, myIngredients, 1, hasThisType<1>()) ); 
		
		// //initialize and run
		// taskmanager2.initialize();
		// taskmanager2.run(1);
		// // taskmanager2.cleanup();
        
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
