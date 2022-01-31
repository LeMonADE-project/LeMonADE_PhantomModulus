/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by
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

#ifndef LEMONADE_UPDATER_SETUP_TMDOUBLESTARSAB
#define LEMONADE_UPDATER_SETUP_TMDOUBLESTARSAB
/**
 * @file
 *
 * @class UpdaterAddTMDoubleStars
 *
 * @brief Updater to create a solution of monodisperse branched stars.
 *
 * @details This is a simple implementation of a system setup starting from an empty ingredients
 * or a system with some monomers inside. This updater requires FeatureAttributes.
 * Two tags are added to the monomers in alternating manner, usually needed for GPU computing.  // not realized in this first attempt
 *                                                                                              // RS, 05 Nov 2019, but presumably 29 July 2020
 *
 * @tparam IngredientsType
 *
 * @param ingredients_ The system, holding either an empty simulation box for system setup
 * or a prefilled ingredients where the stars shall be added
 * @param NStar_ number of stars that are added to ingredients
 * @param NMonoPerStar_ number of monomers in each star
 * @param NMonoPerBranch_ number of monomers in each branch (excluding central monomer)
 * @param NBranchPerStar_ number of branches in each star
 * @param type1_ attribute tag of "even" monomers
 * @param type2_ attribute tag of "odd" monomers
 **/

// #include <LeMonADE/updater/UpdaterAbstractCreate.h>
#include <LeMonADE/utility/Vector3D.h>
#include <cmath>

#include <LeMonADE_PM/updater/UpdaterAbstractCreateAllBondVectors.h>

template<class IngredientsType>
class UpdaterAddTMDoubleStars: public UpdaterAbstractCreateAllBonds<IngredientsType>
{
  typedef UpdaterAbstractCreateAllBonds<IngredientsType> BaseClass;

public:
  UpdaterAddTMDoubleStars(IngredientsType& ingredients_, uint32_t NStar_, 
                  uint32_t NMonoPerBranch_, uint32_t nRings_, uint32_t NBranchPerStar_);

  virtual void initialize();
  virtual bool execute();
  virtual void cleanup();

  //! getter function for number of stars
  const int32_t getNStar() const {return NStar;}

  //! getter function for number of monomers in branch
  const int32_t getNMonoPerBranch() const {return NMonoPerBranch;}

  //! getter function for number of branches in stars
  const int32_t getNBranchPerStar() const {return NBranchPerStar;}

  //! getter function for calculated density
  const double getDensity() const {return density;}

private:
  // provide access to functions of UpdaterAbstractCreate used in this updater
  using BaseClass::ingredients;
  using BaseClass::addMonomerToParent;
  using BaseClass::addSingleMonomer;
  using BaseClass::linearizeSystem;
  double prob_q(double n, double N, double m) {
        return (m+1.)/(N-n-1.);
    }
    //! number of rings 
    uint32_t nRings;
    //! number of stars in the box
    uint32_t NStar;

    //! number of monomers in a branch
    uint32_t NMonoPerBranch;  

    //! number of branches in a star
    uint32_t NBranchPerStar;

    //! lattice occupation density
    double density;

    //! bool for execution
    bool wasExecuted;

    // adds a chain of with nMonomers monomers to the parentID 
    void createChain(uint parentID, uint32_t nMonomers);
    // add a chain to the system at a random positions 
    void createChain(uint32_t nMonomers);

    //! 
    const uint32_t steps;
    //! 
    std::vector<uint32_t> invCPF;
    // random number generator (globally seeded)
    RandomNumberGenerators rng;
};

/**
* @brief Constructor handling the new systems paramters
* 
* @param ingredients_ a reference to the IngredientsType - mainly the system
* @param NStar_ number of double stars to be added in the system instead of solvent
* @param NBranchPerStar_ number of branches in each star
* @param NMonoPerBranch_ number of monomers in each branch (excluding central monomer)
*/


template < class IngredientsType >
UpdaterAddTMDoubleStars<IngredientsType>::UpdaterAddTMDoubleStars(
    IngredientsType& ingredients_, 
    uint32_t NStar_, 
    uint32_t NMonoPerBranch_, 
    uint32_t nRings_,
    uint32_t NBranchPerStar_
    ):
    BaseClass(ingredients_), NStar(NStar_), NMonoPerBranch(NMonoPerBranch_), nRings(nRings_), NBranchPerStar(NBranchPerStar_), 
    density(0.0), wasExecuted(false),steps(50000)
    {
        invCPF.resize(steps-1,0);
    }

/**
* @brief initialise function, calculate the target density to compare with at the end.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAddTMDoubleStars<IngredientsType>::initialize(){
  std::cout << "initialize UpdaterAddTMDoubleStars" << std::endl;


		//probability distribution function 
        std::cout << "Calculate PDF" << std::endl;
		std::vector<double>  PF((NMonoPerBranch-nRings),0);
        double  sum(0.);
		for (auto i = 1; i < (NMonoPerBranch- nRings); i++ ){
            PF[i]=prob_q(i,NMonoPerBranch,nRings)*(1.-sum);
            // std::cout << "probabilityDisti "<<i<< " "<< PF[i]<< std::endl; 
            sum+=PF[i];
		}
        //calculate convolution of the probability distributions 
        //convolution distribution function 
        std::cout << "Calculate convoluted PDF" << std::endl;
		std::vector<double>  convPF((NMonoPerBranch- nRings)*2,0);
        for (auto i=0; i < 2*(NMonoPerBranch-nRings);i++){
            for(auto j=0; j < i ; j++){
                convPF[i]+=PF[j]*PF[i-j];
            }
            // std::cout << "ConvprobabilityDisti "<<i<< " "<< convPF[i]<< std::endl;
        }
        //calculate the cummulative distribution function 
        //cummulative distribution function 
        std::cout << "Calculate cummulative convoluted PDF" << std::endl;
		std::vector<double> CPF((NMonoPerBranch- nRings)*2,0);
		 for (auto i=1; i < 2*(NMonoPerBranch-nRings);i++){
            CPF[i]+=CPF[i-1]+convPF[i];
        }
        //calculate the inverse cummulative convoluted probability distribution 
        std::cout << "Calculate inverse cummulative convoluted PDF" << std::endl;
		for (auto i =1 ; i < steps-1; i++) {
			uint32_t n=0 ;
			auto prob=static_cast<double> (i) / static_cast<double>(steps ); 
			while (  prob > CPF [ n ]  ){ n++;}
			invCPF[i]=(n-1) + static_cast<uint32_t>(round( (n-(n-1)) *(  prob- CPF[n-1] )/(CPF[n]-CPF[n-1]) ));
			// std::cout << "invCPF: " <<  prob << " " << invCPF[i] <<" "<< n-1 <<" "<<CPF [ n-1 ]<<" "<< n <<" "<<CPF [ n ]<< std::endl;
		}

  execute();
}


template < class IngredientsType >
void UpdaterAddTMDoubleStars<IngredientsType>::createChain(uint32_t parentID, uint32_t nMonomers){
    auto ID(parentID);
    for (auto i=0; i< nMonomers ; i++){
        addMonomerToParent(parentID,1);
        parentID=(ingredients.getMolecules().size()-1);
    }
}

template < class IngredientsType >
void UpdaterAddTMDoubleStars<IngredientsType>::createChain(uint32_t nMonomers){
    addSingleMonomer(1);
    createChain(ingredients.getMolecules().size()-1,nMonomers-1);
}
/**
* @brief Execution of the system creation
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
bool UpdaterAddTMDoubleStars<IngredientsType>::execute(){
    if(wasExecuted) return true;
    std::cout << "execute UpdaterAddTMDoubleStars" << std::endl;
    //initial number of monomer in the system 
    auto nMonomers(ingredients.getMolecules().size());
    //keeps track on the monomers which are added during the procedure 
    auto nAddedMonomers(0);
    for(uint32_t i=0;i<NStar;i++){
        auto mainChainLength(invCPF[rng.r250_rand32() % steps ]);
        nAddedMonomers+=mainChainLength;
        std::cout << "UpdaterAddTMDoubleStars: mainChainLength " << mainChainLength <<std::endl;
        //the next two variables are needed for adding NBranchPerStar-1 chains to this monomers 
        //chain start 
        auto start(0);
        if (ingredients.getMolecules().size()>0 )
            start=(ingredients.getMolecules().size()-1);
        std::cout << "Start: " <<  start <<std::endl; 
        createChain(mainChainLength);
        //chain end
        auto end(ingredients.getMolecules().size()-1);
        std::cout << "End: " <<  end <<std::endl; 
        for (uint32_t j=0; j < NBranchPerStar-1; j++){
            auto chainLength(invCPF[rng.r250_rand32() % steps ]+1);
            nAddedMonomers+=chainLength;
            createChain(start,chainLength);
            ingredients.modifyMolecules()[ ingredients.getMolecules().size()-1 ].setMovableTag(false);  
        }
        for (uint32_t j=0; j < NBranchPerStar-1; j++){
            auto chainLength(invCPF[rng.r250_rand32() % steps ]+1);
            nAddedMonomers+=chainLength;
            createChain(end,chainLength);
            ingredients.modifyMolecules()[ ingredients.getMolecules().size()-1 ].setMovableTag(false);  
        }
    }

    ingredients.synchronize();
    if(nAddedMonomers != (ingredients.getMolecules().size()-nMonomers)){
        std::cout << "initial number of monomers: "<< nMonomers << std::endl;
        std::cout << "final number of monomers: "<< ingredients.getMolecules().size() << std::endl;
        std::cout << "difference should be: "<< nAddedMonomers << " but is "<<  (ingredients.getMolecules().size()-nMonomers)<< std::endl;
        throw std::runtime_error("UpdaterAddTMDoubleStars: number of monomers in molecules does not match the calculated number of monomers!");
    }else{
        std::cout << "number of monomers in the system  =" << (ingredients.getMolecules().size())<<std::endl;
        std::cout << "added number of monomers to the system  =" << (ingredients.getMolecules().size()-nMonomers)<<std::endl;
        wasExecuted=true;
        return true;
    }
}

/**
* @brief Standard clean up.
*
* @tparam IngredientsType Features used in the system. See Ingredients.
*/
template < class IngredientsType >
void UpdaterAddTMDoubleStars<IngredientsType>::cleanup(){

}


#endif /* LEMONADE_UPDATER_SETUP_TMDOUBLESTARSAB */
