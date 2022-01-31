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
#ifndef LEMONADE_PM_UPDATER_UPDATERREADCROSSLINKCONNECTIONSTENDOMER_H 
#define LEMONADE_PM_UPDATER_UPDATERREADCROSSLINKCONNECTIONSTENDOMER_H
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <math.h>
#include <algorithm>
#include <LeMonADE/updater/AbstractUpdater.h>
#include <LeMonADE/utility/DistanceCalculation.h>

/**
 * @class UpdaterReadCrosslinkConnectionsTendomer
 * @brief reads in connections for a system made of linear chains and crosslinker 
 * 
 * @details 
 *   -the input file needs to have a format like:
 *     #Time, ChainID, MonID1, P1X, P1Y, P1Z, MonID2, P2X, P2Y, P2Z
 *   -the chains are before the crosslinks in the bfm file
 *   -the chain length must be at least 1 
 * 
 */

template <class IngredientsType>
class UpdaterReadCrosslinkConnectionsTendomer : public AbstractUpdater
{
public:
    UpdaterReadCrosslinkConnectionsTendomer(
        IngredientsType& ing_, 
        const std::string input_, 
        const double stepwidth_, 
        const double minConversion_): 
        ing(ing_), 
        input(input_), 
        stepwidth(stepwidth_), 
        minConversion(minConversion_),
        nExecutions(0){};
    virtual void initialize();
    virtual bool execute();
    virtual void cleanup(){};

private:
  //! container storing system information about monomers
  IngredientsType& ing;

  //! container storing system information about monomers
  IngredientsType initialIng;

  //! input of the cross link positions and so on..
  const std::string input;

  //!conversion step used in the execute function to increase the conversion 
  const double stepwidth;

  //! minimal conversion to be read in before further analyse
  double minConversion;

  //!number of maximum connections for a cross link
  uint32_t NMaxConnection;

  //!numbe of monomers per chain 
  uint32_t NMonomerPerChain;

  //!number of executions;
  uint32_t nExecutions;

  //! connects the cross link to the chain
  bool ConnectCrossLinkToChain(uint32_t MonID, uint32_t chainID);
  
  //!bond Table 
  std::map<std::pair<uint32_t,uint32_t>,std::vector<uint32_t> > bondTable;
};
template <class IngredientsType>
bool  UpdaterReadCrosslinkConnectionsTendomer<IngredientsType>::ConnectCrossLinkToChain(uint32_t MonID, uint32_t chainID){
    std::pair<uint32_t,uint32_t> key(MonID,chainID-1);
    if (bondTable.find(key) != bondTable.end()){
        if ( !ing.getMolecules().areConnected( MonID,bondTable.at(key)[0] ) ){
            ing.modifyMolecules().connect(MonID,bondTable.at(key)[0] );
            return true ; 
        }
        //chain can only have one neighbor (without this statement a more or less random partner would be connected to the structure !!)
        if ( bondTable.at(key).size()>1 ) {
            ing.modifyMolecules().connect(MonID,bondTable.at(key)[1] );
            return true ; 
        }
    }
    return false; 
}
/**
 * @brief read in connections up to the minimum conversion given 
 * */
template <class IngredientsType>
void UpdaterReadCrosslinkConnectionsTendomer<IngredientsType>::initialize(){
    //assume a stochiometric mixture
    NMaxConnection=4*ing.getNumCrossLinkers();
    NMonomerPerChain = 2*ing.getNumMonomersPerChain();
    std::cout << "Number of maximum connection: " << NMaxConnection << std::endl;
    //erase bonds between reactive monomers
    for (uint32_t i =0 ; i <  ing.getMolecules().size(); i++)
        if (ing.getMolecules()[i].isReactive() )
            for(size_t j=0; j < ing.getMolecules().getNumLinks(i);j++){
                uint32_t neighbor(ing.getMolecules().getNeighborIdx(i,j));
                if (ing.getMolecules()[neighbor].isReactive()){
                    ing.modifyMolecules().disconnect(i, neighbor );
                    uint32_t chainMonomer(std::min(i,neighbor) );
                    uint32_t crosslink(std::max(i,neighbor));
                    uint32_t chainID( (chainMonomer-chainMonomer%NMonomerPerChain)/NMonomerPerChain);
                    bondTable[std::pair<uint32_t,uint32_t>(crosslink,chainID) ].push_back(chainMonomer) ;
                }
            }
    std::cout << "Erase " << bondTable.size() << " bonds." <<std::endl;
    ing.synchronize();
    initialIng=ing;
}
/**
 * @brief read in connections up tp the next step 
 * @details Copies the initial ingredients to the current one and adds connections up to the current conversion. 
 * */
template <class IngredientsType>
bool UpdaterReadCrosslinkConnectionsTendomer<IngredientsType>::execute(){
    //reset the ingredients container to the inital one
    ing = initialIng;
    //stream reading input file
    std::ifstream stream;
    stream.open(input);
    if (stream.fail())
      throw std::runtime_error(std::string("error opening input file ") + input + std::string("\n"));
    //current conversion 
    auto conversion = minConversion + static_cast<double>(nExecutions) * stepwidth;
    std::cout << "Current conversion is " <<conversion <<std::endl;
    //read in number of lines to reach the current conversion 
    //(always start from the initial state and thus read in connectiosn from the beginning)
    uint32_t ReadNLines = (floor(NMaxConnection * conversion));
    //counter for the new connections
    uint32_t NewConnections(0);
    std::cout << "Start reading " << ReadNLines << " number of lines " << std::endl;
    while (NewConnections < ReadNLines && stream.good()){
        std::string line;
        getline(stream, line);
        if (line.empty() || line.at(0) == '#' )
            continue;
        std::stringstream ss;
        uint32_t Time, createBreak, ChainID, nSegments, MonID1, P1X, P1Y, P1Z, MonID2, P2X, P2Y, P2Z;
        ss << line;
        ss >> Time >> createBreak>>  ChainID >> nSegments >>MonID1 >> P1X >> P1Y >> P1Z >> MonID2 >> P2X >> P2Y >> P2Z;
        #ifdef DEBUG
            std::cout << ss.str() << std::endl;
        #endif //DEBUG//
        ing.modifyMolecules().setAge(Time);
        if ( ConnectCrossLinkToChain(MonID1, ChainID+1) ) {
        }else if( ConnectCrossLinkToChain(MonID2, ChainID+1) ) {
        }else {
            std::stringstream errormessage;
            errormessage << "There was no such a connection in the bfm file for monomer ID= " << MonID1  <<" of with ID=" << MonID2 <<  " with chainID=" << ChainID<< "\n";
            throw std::runtime_error(errormessage.str());
        }
        NewConnections++;
    }
    std::cout << "Read and add " << NewConnections << "/" << NMaxConnection 
              << " to the system at time "
              << ing.getMolecules().getAge() << std::endl;
    ing.synchronize();
    nExecutions++;
    std::cout << "UpdaterReadCrosslinkConnectionsTendomer::execute " << nExecutions << " times.\n";
    //close the filestream and return false if the file has ended and thus the updater has nothing more to do
    if (stream.eof()) {
        stream.close();
        return false;
    }else{
        stream.close();
        return true;
    }
}

#endif