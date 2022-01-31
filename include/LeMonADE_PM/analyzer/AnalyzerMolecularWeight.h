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

#ifndef LEMONADE_PM_ANALYZER_MOLECULAR_WEIGTH_H
#define LEMONADE_PM_ANALYZER_MOLECULAR_WEIGTH_H

#include <string>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>
#include <LeMonADE/utility/MonomerGroup.h>
#include <LeMonADE/utility/DepthIterator.h>
#include <LeMonADE/analyzer/AnalyzerAbstractDump.h>


/*************************************************************************
 * definition of AnalyzerMolecularWeight class
 * ***********************************************************************/

/**
 * @file
 *
 * @class AnalyzerMolecularWeight
 *
 * @tparam IngredientsType Ingredients class storing all system information( e.g. monomers, bonds, etc).
 * 
 * @details Calculates the number and weight averaged molecular weight.
 *
 */
template < class IngredientsType > 
class AnalyzerMolecularWeight : public AnalyzerAbstractDump<IngredientsType,double>
{
    typedef AnalyzerAbstractDump<IngredientsType,double> BaseClass; 
private:
	//! typedef for the underlying container holding the monomers
	typedef typename IngredientsType::molecules_type molecules_type;

    std::string basename;
    
protected:
  using BaseClass::ingredients;
	using BaseClass::Data;
	using BaseClass::MCSTimes;
	using BaseClass::dumpTimeSeries;
    
public:
	//! constructor
	AnalyzerMolecularWeight(const IngredientsType& ing, std::string fileSuffix_);
	//! destructor. does nothing
	virtual ~AnalyzerMolecularWeight(){}
	//! calculate the moleculare weight distribution
	virtual void initialize();
  virtual bool execute();
};

/*************************************************************************
 * implementation of memebers
 * ***********************************************************************/

/**
 * @param ing reference to the object holding all information of the system
 * @param fileSuffix output file name. defaults to "Rg2TimeSeries.dat".
 * */
template<class IngredientsType>
AnalyzerMolecularWeight<IngredientsType>::AnalyzerMolecularWeight(
	const IngredientsType& ing,std::string fileSuffix_)
:BaseClass(ing,fileSuffix_) {}

template <class IngredientsType>
void AnalyzerMolecularWeight<IngredientsType>::initialize()
{
  basename=BaseClass::getOutputFilename() ; 
  BaseClass::setBufferSize(1);
  BaseClass::setNumberOfColumns(1);
  execute();
}
template< class IngredientsType >
bool AnalyzerMolecularWeight<IngredientsType>::execute()
{
    double conversion = ingredients.getConversion();
    std::stringstream comment;
    comment << "Created by AnalyzerMolecularWeight\n"
            << "conversion=" << conversion << "%";
    BaseClass::setComment(comment.str());

    //calculate the molecular weight distribution 
    std::vector<MonomerGroup<molecules_type> > groups;
    fill_connected_groups(ingredients.getMolecules(), groups, MonomerGroup<molecules_type>(ingredients.getMolecules()), alwaysTrue());
    //key is the size of the molecule and the value is the number of occurence 
    std::map<uint32_t,uint32_t> dist;
    for (auto it=groups.begin(); it!=groups.end();it++ )
      dist[it->size()]++;
    for ( auto&  it : dist )
    {
      Data[0].push_back(it.second);
      MCSTimes.push_back(it.first);
    }
    std::stringstream output; 
    output << "MolecularWeight_c" << conversion <<  basename<< ".dat";
    BaseClass::resetIsFirstFileDump();
    BaseClass::setOutputFilename(output.str());
    dumpTimeSeries();
}
#endif
