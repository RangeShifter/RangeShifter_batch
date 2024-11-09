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

 There is ONE instance of a Population for each Species within each Patch
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
	Population(); // default constructor
	// constructor for a Population of a specified size
	Population(Species* pSpecies, Patch* pPatch, int ninds, int	resol);
	~Population();

	traitsums getIndTraitsSums(Species*);
	popStats getStats();
	Patch* getPatch() { return pPatch; }
	Species* getSpecies();
	int getNInds();
	int totalPop();
	// return no. of Individuals in a specified stage
	int stagePop(int stg);
	void localExtinction(int option);
		// option: 	0 - random local extinction probability
		//			1 - local extinction probability gradient

	void extirpate(); // Remove all individuals
	void reproduction(
		const float,	// local carrying capacity
		const float,	// effect of environmental gradient and/or stochasticty
		const int			// Landscape resolution
	);
	// Following reproduction of ALL species, add juveniles to the population
	void fledge();

	// Determine which individuals will disperse
	void emigration(float localK);
	// All individuals emigrate after patch destruction
	void allEmigrate();

	// If an individual has been identified as an emigrant, remove it from the Population
	disperser extractDisperser(int ix);
	// For an individual identified as being in the matrix population:
	// if it is a settler, return its new location and remove it from the current population
	// otherwise, leave it in the matrix population for possible reporting before deletion
	disperser extractSettler(int ix);

	// Add a specified individual to the population
	void recruit(Individual* pInd);
	Individual* sampleInd() const;
	void sampleIndsWithoutReplacement(string n, const set<int>& sampleStages);
	int sampleSize() const;
	vector<Individual*> getIndividualsInStage(int stage);

	int transfer( // Executed for the Population(s) in the matrix only
		Landscape* pLandscape,
		short landIx,
		short nextseason
	);

	// Determine whether there is a potential mate present in a patch which a potential
	// settler has reached
	bool isMatePresent(Cell* pCell, short sex);

	// Determine survival and development and record in individual's status code
	// Changes are NOT applied to the Population at this stage
	void drawSurvivalDevlpt(
		bool resolveJuvs, 
		bool resolveAdults,
		bool resolveDev, 
		bool resolveSurv
	);
	void applySurvivalDevlpt(); // Apply survival changes to the population
	void ageIncrement();

	// Open population file and write header record
	bool outPopHeaders(int landNr, bool patchModel);

	void outPopulation( // Write record to population file
		int rep,
		int yr,
		int gen,
		bool envLocal,
		float eps,
		bool patchModel,
		bool writeEnv,
		bool gradK
	);

	// Open individuals file and write header record
	void outIndsHeaders(int rep, int landnr, bool patchModel);
	// Write records to individuals file
	void outIndividual(Landscape* pLandscape, int rep, int yr, int gen);
	void outputTraitPatchInfo(ofstream& outtraits, int rep, int yr, int gen, bool patchModel);
	traitsums outTraits(ofstream& outtraits, const bool& writefile);
	void outputGeneValues(ofstream& ofsGenes, const int& yr, const int& gen) const;
	
	void clean(); // Remove zero pointers to dead or dispersed individuals

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
	Species* pSpecies;
	Patch* pPatch;
	int nInds[gMaxNbStages][gMaxNbSexes];

	vector<Individual*> inds; // all individuals in population except ...
	vector<unique_ptr<Individual>> juvs; // ... juveniles until reproduction of ALL species
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

