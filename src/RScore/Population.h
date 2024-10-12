/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2020 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Damaris Zurell
 *
 *	This file is part of RangeShifter.
 *
 *	RangeShifter is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	RangeShifter is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with RangeShifter. If not, see <https://www.gnu.org/licenses/>.
 *
 --------------------------------------------------------------------------*/


 /*------------------------------------------------------------------------------

 RangeShifter v2.0 Population

 Implements the Population class

 There is ONE instance of a Population for each Species within each SubCommunity
 (including the matrix). The Population holds a list of all the Individuals in
 the Population.

 The matrix Population(s) hold(s) Individuals which are currently in the process
 of transfer through the matrix.

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 25 June 2021 by Steve Palmer

 ------------------------------------------------------------------------------*/

#ifndef PopulationH
#define PopulationH

#include <vector>
#include <algorithm>
using namespace std;

#include "Parameters.h"
#include "Individual.h"
#include "Species.h"
#include "Landscape.h"
#include "Patch.h"
#include "Cell.h"
#include "NeutralStatsManager.h"
#include "Population.h"

class Patch;

//---------------------------------------------------------------------------

struct popStats {
	Species* pSpecies; 
	Patch* pPatch; 
	int spNum, nInds, nNonJuvs, nAdults; 
	bool breeding;
};
struct disperser {
	Individual* pInd; 
	Cell* pCell; 
	bool yes;
};

struct traitsums {
	vector<int> ninds = vector<int>(gMaxNbSexes, 0);
	vector<double> sumD0 = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqD0 = vector<double>(gMaxNbSexes, 0);
	vector<double> sumAlpha = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqAlpha = vector<double>(gMaxNbSexes, 0);
	vector<double> sumBeta = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqBeta = vector<double>(gMaxNbSexes, 0);
	vector<double> sumDist1 = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqDist1 = vector<double>(gMaxNbSexes, 0);
	vector<double> sumDist2 = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqDist2 = vector<double>(gMaxNbSexes, 0);
	vector<double> sumProp1 = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqProp1 = vector<double>(gMaxNbSexes, 0);
	vector<double> sumDP = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqDP = vector<double>(gMaxNbSexes, 0);
	vector<double> sumGB = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqGB = vector<double>(gMaxNbSexes, 0);
	vector<double> sumAlphaDB = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqAlphaDB = vector<double>(gMaxNbSexes, 0);
	vector<double> sumBetaDB = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqBetaDB = vector<double>(gMaxNbSexes, 0);
	vector<double> sumStepL = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqStepL = vector<double>(gMaxNbSexes, 0);
	vector<double> sumRho = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqRho = vector<double>(gMaxNbSexes, 0);
	vector<double> sumS0 = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqS0 = vector<double>(gMaxNbSexes, 0);
	vector<double> sumAlphaS = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqAlphaS = vector<double>(gMaxNbSexes, 0);
	vector<double> sumBetaS = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqBetaS = vector<double>(gMaxNbSexes, 0);
	vector<double> sumGeneticFitness = vector<double>(gMaxNbSexes, 0);
	vector<double> ssqGeneticFitness = vector<double>(gMaxNbSexes, 0);
};

class Population {

public:
	Population(void); // default constructor
	Population( // constructor for a Population of a specified size
		Species*,	// pointer to Species
		Patch*,		// pointer to Patch
		int,			// no. of Individuals
		int				// Landscape resolution
	);
	~Population(void);
	traitsums getIndTraitsSums(Species*);
	popStats getStats(void);
	Species* getSpecies(void);
	int getNInds(void);
	int totalPop(void);
	int stagePop( // return no. of Individuals in a specified stage
		int	// stage
	);
	void extirpate(void); // Remove all individuals
	void reproduction(
		const float,	// local carrying capacity
		const float,	// effect of environmental gradient and/or stochasticty
		const int			// Landscape resolution
	);
	// Following reproduction of ALL species, add juveniles to the population
	void fledge(void);
	void emigration( // Determine which individuals will disperse
		float   // local carrying capacity
	);
	void allEmigrate(void); // All individuals emigrate after patch destruction
	// If an individual has been identified as an emigrant, remove it from the Population
	disperser extractDisperser(
		int		// index no. to the Individual in the inds vector
	);
	// For an individual identified as being in the matrix population:
	// if it is a settler, return its new location and remove it from the current population
	// otherwise, leave it in the matrix population for possible reporting before deletion
	disperser extractSettler(
		int   // index no. to the Individual in the inds vector
	);
	void recruit( // Add a specified individual to the population
		Individual*	// pointer to Individual
	);
	Individual* sampleInd() const;
	void sampleIndsWithoutReplacement(string n, const set<int>& sampleStages);
	int sampleSize() const;
	vector<Individual*> getIndividualsInStage(int stage);
#if RS_RCPP
	int transfer( // Executed for the Population(s) in the matrix only
		Landscape*,	// pointer to Landscape
		short,				// landscape change index
		short				// year
	);
	// Determine whether there is a potential mate present in a patch which a potential
	// settler has reached
	bool matePresent(
		Cell*,	// pointer to the Cell which the potential settler has reached
		short		// sex of the required mate (0 = female, 1 = male)
	);
#else
	int transfer( // Executed for the Population(s) in the matrix only
		Landscape*,	// pointer to Landscape
		short				// landscape change index
	);
	// Determine whether there is a potential mate present in a patch which a potential
	// settler has reached
	bool matePresent(
		Cell*,	// pointer to the Cell which the potential settler has reached
		short		// sex of the required mate (0 = female, 1 = male)
	);
#endif // RS_RCPP
	// Determine survival and development and record in individual's status code
	// Changes are NOT applied to the Population at this stage
	void survival0(
		float,	// local carrying capacity
		short,	// option0:	0 - stage 0 (juveniles) only
		//	  			1 - all stages
		//					2 - stage 1 and above (all non-juveniles)
		short 	// option1:	0 - development only (when survival is annual)
						//	  	 		1 - development and survival
						//	  	 		2 - survival only (when survival is annual)
	);
	void survival1(void); // Apply survival changes to the population
	void ageIncrement(void);
	bool outPopHeaders( // Open population file and write header record
		int,	// Landscape number (-999 to close the file)
		bool	// TRUE for a patch-based model, FALSE for a cell-based model
	);
	void outPopulation( // Write record to population file
		int,		// replicate
		int,		// year
		int,		// generation
		float,	// epsilon - global stochasticity value
		bool,		// TRUE for a patch-based model, FALSE for a cell-based model
		bool,		// TRUE to write environmental data
		bool		// TRUE if there is a gradient in carrying capacity
	);

	void outIndsHeaders( // Open individuals file and write header record
		int,	// replicate
		int,	// Landscape number (-999 to close the file)
		bool	// TRUE for a patch-based model, FALSE for a cell-based model
	);
	void outIndividual( // Write records to individuals file
		Landscape*,	// pointer to Landscape
		int,				// replicate
		int,				// year
		int,				// generation
		int					// Patch number
	);

	void outputTraitPatchInfo(ofstream& outtraits, int rep, int yr, int gen, bool patchModel);
	traitsums outTraits(ofstream& outtraits);

	void outputGeneValues(ofstream& ofsGenes, const int& yr, const int& gen) const;
	void clean(void); // Remove zero pointers to dead or dispersed individuals

	void updatePopNeutralTables();
	double getAlleleFrequency(int locus, int allele);
	int getAlleleTally(int locus, int allele);
	int getHeteroTally(int locus, int allele);
	int countHeterozygoteLoci();
	vector<int> countNbHeterozygotesEachLocus();
	double computeHs();

#ifndef NDEBUG
	// Testing only
	void clearInds() { inds.clear(); } // empty inds vector to avoid deallocating inds when used in test
#endif // NDEBUG

private:
	short nStages;
	short nSexes;
	Species* pSpecies;	// pointer to the species
	Patch* pPatch;			// pointer to the patch
	int nInds[gMaxNbStages][gMaxNbSexes];		// no. of individuals in each stage/sex

	vector<Individual*> inds; // all individuals in population except ...
	vector<Individual*> juvs; // ... juveniles until reproduction of ALL species
	// has been completed

	vector<Individual*> sampledInds;
	vector<NeutralCountsTable> popNeutralCountTables;
	void resetPopNeutralTables();
};

//---------------------------------------------------------------------------

extern paramGrad* paramsGrad;
extern paramStoch* paramsStoch;
extern paramInit* paramsInit;
extern paramSim* paramsSim;
extern RSrandom* pRandom;

//---------------------------------------------------------------------------
#endif

