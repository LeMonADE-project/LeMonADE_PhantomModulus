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

#ifndef LEMONADE_PM_UPDATER_MOVES_MOVENONLINEARFORCEEQUILIBRIUM_H
#define LEMONADE_PM_UPDATER_MOVES_MOVENONLINEARFORCEEQUILIBRIUM_H
#include <limits>
#include <fstream>
#include <sys/stat.h>
#include <LeMonADE_PM/updater/moves/MoveForceEquilibriumBase.h>
#include <LeMonADE/utility/DistanceCalculation.h>
#include <LeMonADE_PM/utility/neighborX.h>

/*****************************************************************************/
/**
 * @file
 *
 * @class MoveNonLinearForceEquilibrium
 *
 * @brief Standard local bfm-move on simple cubic lattice for the scBFM.
 *
 * @details The class is a specialization of MoveLocalBase using the (CRTP) to avoid virtual functions.
 * Here implemented for networks made of monodisperse chains.
 **/
/*****************************************************************************/

class MoveNonLinearForceEquilibrium:public MoveForceEquilibriumBase<MoveNonLinearForceEquilibrium>{
public:
    //! constructor for the class MoveNonLinearForceEquilibrium taking the filename of the force extension relation and the realaxation parameter for the chain 
    MoveNonLinearForceEquilibrium(std::string filename_="", double relaxationChain_=1.):
        filename(filename_),bondlength(2.68), accuracy(.1) {
            if( !filename.empty() ) createTable();
            setRelaxationParameter(relaxationChain_);
        };
    // overload initialise function to be able to set the moves index and direction if neccessary
    template <class IngredientsType> void init(const IngredientsType& ing);
    template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index);
    template <class IngredientsType> void init(const IngredientsType& ing, uint32_t index, VectorDouble3 dir );

    template <class IngredientsType> bool check(IngredientsType& ing);
    template< class IngredientsType> void apply(IngredientsType& ing);

    //!read file to create a lookup table of the force extension curve 
    void createTable();
    // //!read file to create a lookup table of the force extension curve 
    // void createTable(std::string filename_){
    //     setFilename(filename); 
    //     createTable();}
    //! set the filename for the force extension data 
    void setFilename(std::string filename_){filename=filename_;createTable();}
    //! get the filename for the force extension data 
    std::string const getFilename(){return filename;}
    
    //! set the relaxation parameter for the cross link
    void setRelaxationParameter(double relaxationChain_){
        relaxationChain=relaxationChain_;
        springConstant = relaxationChain*bondlength*bondlength/(3.);
        std::cout << "Set relaxation parameter to " << relaxationChain << " corresponding to a spring constant " <<springConstant <<std::endl;
    }
    //! get the relaxation parameter for the cross link 
    double getRelaxationParameter(){return relaxationChain;} 

    //uses the read force extension relation 
    VectorDouble3 EF(VectorDouble3 extensionVector){
        if (extensionVector == VectorDouble3(0.,0.,0.) )
            return VectorDouble3(0.,0.,0.);
        double length( extensionVector.getLength() );

        #ifdef DEBUG
        if ( max_extension < length  ){
            std::stringstream errormessage; 
            errormessage << "The length of the extension vector is greater than the maximum given in the file: \n"
                         << "length is  " << length 
                         << " maximum is " << max_extension << "\n";
            throw std::runtime_error(errormessage.str());
        }
        if (length == 0 )
            return VectorDouble3(0.,0.,0.);
        #endif
        #ifdef DEBUG
        std::cout <<length << "\t" << force_extension[static_cast<uint32_t>(round(length/accuracy))] << "\t" << EFGauss(extensionVector).getLength() << "\t" 
         << extensionVector << "\t" 
         <<  EFGauss(extensionVector) << "\t"<< extensionVector.normalize()*(force_extension[static_cast<uint32_t>(round(length/accuracy))])<<"\t"  
        << "\n"; 
        #endif 
        if (length == 0)
            return VectorDouble3(0.,0.,0.);
        auto x(length/accuracy);
        auto up (static_cast<uint32_t>(ceil(x)+1 ));
        auto down(static_cast<uint32_t>(floor(x)));
        auto amplitude=force_extension[down] + (force_extension[up]-force_extension[down])/static_cast<double>(up-down)*(x-down)/static_cast<double>(up-down);
        return extensionVector.normalize()*(amplitude);
    }
    //Gaussina force extension relation 
    //f=R*3/(N*b^2)
    VectorDouble3 EFGauss(VectorDouble3 extensionVector){
        return extensionVector/springConstant;
    }
    //Gaussian extension force relation 
    //R=-f/3*N*b^2
    VectorDouble3 FE(VectorDouble3 const force ) const {
        return force*springConstant;
    }

private:
    //average square bond length 
    const double bondlength;
    //filename filename for the force extension curve 
    std::string filename;
    
    //minimum force in the file 
    double min_force;
    //maximum force in the file 
    double max_force;
    //force steps 
    double force_steps;
   
    //min_extension in the file 
    double min_extension;
    //maximum extensions in the file 
    double max_extension;
    //extension steps 
    double extension_steps;

    //the accuracy is the steps in which the extensions vector grows
    const double accuracy; 

    //spring constant for the equivalent chain used for the relaxation of the cross links 
    double springConstant;
   
    //force extension mapping       
    std::map<double, double> extension_force;
    //extension force mapping index=extension rounded to int 
    std::vector<double> force_extension;
    //an equivalent chain which relaxes the cross link
    double relaxationChain;
    //calculate the shift for the cross link
    template< class IngredientsType >
    VectorDouble3 CalculateShift(IngredientsType& ing ){
        std::vector<neighborX> Neighbors(ing.getCrossLinkNeighborIDs(this->getIndex()) );
        VectorDouble3 force(0.,0.,0.);
        VectorDouble3 shift(0.,0.,0.);
        double avNSegments(0.);
        int32_t number_of_neighbors(Neighbors.size());
        if (number_of_neighbors > 0) {
            VectorDouble3 Position(ing.getMolecules()[this->getIndex()].getVector3D());      
            for (size_t i = 0; i < number_of_neighbors; i++){
                VectorDouble3 vec(ing.getMolecules()[Neighbors[i].ID].getVector3D()-Neighbors[i].jump-Position);
                force+=EF(vec);
            }
            shift=FE(force/(static_cast<double>(number_of_neighbors) ));
        }
        return shift;
    };

    //! check is file exists:
    inline bool fileExists (const std::string& name) {
        struct stat buffer;   
        return (stat (name.c_str(), &buffer) == 0); 
    }
};
/////////////////////////////////////////////////////////////////////////////
/////////// implementation of the members ///////////////////////////////////
void MoveNonLinearForceEquilibrium::createTable(){
    if (fileExists(filename )){
        std::ifstream in(filename);
        uint32_t counter(0);
        extension_force.insert(extension_force.end(),std::pair<double,double>(0., 0.));
        min_force=0.;
        min_extension=0.;
        while(in.good() && in.peek()!=EOF){
            std::string line;
            getline(in, line);
            //ignore comments and blank lines 
            while (line.at(0) == '#' || line.empty() ) //go to next line
                continue;
            //read data 
            double force, extension;
            std::stringstream ss ;
            ss<< line;
            ss>>force >> extension; 
            if(min_force > force ) min_force=force; 
            if(max_force < force ) max_force=force; 
            if(min_extension > extension ) min_extension=extension; 
            if(max_extension < extension ) max_extension=extension; 
            extension_force.insert(extension_force.end(),std::pair<double,double>(force, extension));
            if(counter==1){
                force_steps=max_force-min_force;
            }else if(counter>1)
            counter++;
        }
        in.close();
        //make lookup for the extension force relation 
        //make a entry from 0 to max_extension in steps of 1 
        auto it_last=extension_force.begin();
        for ( auto r=0;r<static_cast<uint32_t>(max_extension/accuracy); r++   ){
            if(r==0)
                force_extension.push_back(0.);
            else{
                for (auto it=it_last; it !=extension_force.end();it++){
                    if(it->second > static_cast<double>(r*accuracy)){
                        //at the force linear interpolated in between the two current forces
                        auto deltaForce(it->first-it_last->first);
                        auto deltaExtension(it->second-it_last->second);
                        auto factor( (static_cast<double>(r*accuracy)-it_last->second)/deltaExtension );
                        //at interpolated force
                        force_extension.push_back(it_last->first+deltaForce*factor);
                        break;
                    }
                    it_last=it;
                }
            }
        }
        std::cout << "MoveNonLinearForceEquilibrium::createTable() force extension" <<std::endl;
        for (auto i=0; i<40; i++ )
            std::cout << "FECurve: " << i << "\t" << i*accuracy<< "\t" << force_extension[i]<<std::endl;

        std::cout << "MoveNonLinearForceEquilibrium::createTable(): \n" 
                << "min_force=" << min_force <<"\n"
                << "max_force=" << max_force <<"\n"
                << "min_extension=" << min_extension <<"\n"
                << "max_extension=" << max_extension <<"\n";
    }else{
        std::cerr<< "Provide a filename for the MoveNonLinearForceEquilibrium!\n" ;
    }
}
/*****************************************************************************/
/**
 * @brief Initialize the move.
 *
 * @details Resets the move probability to unity. Dice a new random direction and
 * Vertex (monomer) index inside the graph.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template <class IngredientsType>
void MoveNonLinearForceEquilibrium::init(const IngredientsType& ing)
{
    this->resetProbability();

    //draw index
    this->setIndex( (this->randomNumbers.r250_rand32()) %(ing.getMolecules().size()) );

    //calculate the shift of the cross link
    this->setShiftVector(CalculateShift(ing));
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be connected
 **/
template <class IngredientsType>
void MoveNonLinearForceEquilibrium::init(const IngredientsType& ing, uint32_t index)
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveNonLinearForceEquilibrium::init(ing, index): index out of range!");

  //calculate the shift of the cross link
  this->setShiftVector(CalculateShift(ing));
  
}

/*****************************************************************************/
/**
 * @brief Initialize the move with a given monomer index.
 *
 * @details Resets the move probability to unity. Dice a new random direction.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @param index index of the monomer to be connected
 * @param bondpartner index of the monomer to connect to 
 **/
template <class IngredientsType>
void MoveNonLinearForceEquilibrium::init(const IngredientsType& ing, uint32_t index, VectorDouble3 dir )
{
  this->resetProbability();

  //set index
  if( (index >= 0) && (index <= (ing.getMolecules().size()-1)) )
    this->setIndex( index );
  else
    throw std::runtime_error("MoveNonLinearForceEquilibrium::init(ing, index, bondpartner): index out of range!");

  //calculate the shift of the cross link
  this->setShiftVector(dir);
}

/*****************************************************************************/
/**
 * @brief Check if the move is accepted by the system.
 *
 * @details This function delegates the checking to the Feature.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 * @return True if move is valid. False, otherwise.
 **/
template <class IngredientsType>
bool MoveNonLinearForceEquilibrium::check(IngredientsType& ing)
{
  //send the move to the Features to be checked
  return ing.checkMove(ing,*this);
}

/*****************************************************************************/
/**
 * @brief Apply the move to the system , e.g. add the displacement to Vertex (monomer) position.
 *
 * @details As first step: all Feature should apply the move using applyMove().\n
 * Second: Modify the positions etc. of the Vertex etc.
 *
 * @param ing A reference to the IngredientsType - mainly the system
 **/
template< class IngredientsType>
void MoveNonLinearForceEquilibrium::apply(IngredientsType& ing)
{

	//move must FIRST be applied to the features
	ing.applyMove(ing,*this);
	//THEN the position can be modified
	ing.modifyMolecules()[this->getIndex()]+=this->getShiftVector();
}

#endif /*MOVEFORCEEQUILIBRIUM_H*/
