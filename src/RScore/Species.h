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

 RangeShifter v2.0 Species

 Implements the Species class

 There is ONE instance of a Species for each species within the Community
 AND THIS IS CURRENTLY LIMITED TO A SINGLE SPECIES.
 The class holds all the demographic and dispersal parameters of the species.

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 28 July 2021 by Greta Bocedi

 ------------------------------------------------------------------------------*/

#ifndef SpeciesH
#define SpeciesH

#include <ranges>
#include <map>
#include <set>
#include <memory>

#include "Parameters.h"
#include "SpeciesTrait.h"
#include "QuantitativeTrait.h"

class SpeciesTrait;

// structures for demographic parameters

struct demogrParams {
	short repType;
	short repSeasons;
	float propMales; 
	float harem; 
	float bc; 
	float lambda;
	bool stageStruct;
};
struct stageParams {
	short nStages;
	short repInterval; 
	short maxAge;
	short survival;
	float probRep;
	bool fecDens;  
	bool fecStageDens; 
	bool devDens; 
	bool devStageDens;
	bool survDens;
	bool survStageDens;
	bool disperseOnLoss;
};
struct densDepParams {
	float devCoeff; 
	float survCoeff;
};

// structures for emigration parameters

struct emigRules {
	bool densDep;
	bool stgDep; 
	bool sexDep; 
	bool indVar;
	short emigStage;
	short emigTrait[2];
};
struct emigTraits {
	float d0; 
	float alpha; 
	float beta;

	emigTraits() : d0(0.0), alpha(0.0), beta(0.0) {};

	emigTraits(const emigTraits& e) : d0(e.d0), alpha(e.alpha), beta(e.beta) {};

	emigTraits* clone() { return new emigTraits(*this); }

	void divideTraitsBy(int i) {

		d0 /= i;
		alpha /= i;
		beta /= i;
	}
};

// structures for transfer parameters

struct transferRules {
	bool usesMovtProc; 
	bool stgDep; 
	bool sexDep;
	bool distMort;
	bool indVar;
	bool twinKern;
	bool habMort;
	short moveType;
	bool usesCosts;
	short movtTrait[2];
};
struct trfrKernelParams {
	float	meanDist1;
	float	meanDist2;
	float	probKern1;
};
struct trfrMortParams {
	float fixedMort;
	float mortAlpha;
	float mortBeta;
};
struct trfrMovtParams {
	short pr;
	short prMethod; 
	short memSize; 
	short goalType;
	float dp; 
	float gb;
	float alphaDB;
	int betaDB;
	float stepMort; 
	float stepLength; 
	float rho;
	bool straightenPath;
};
struct trfrCRWTraits {
	float stepMort;
	float stepLength;
	float rho;
	bool straightenPath;
};
struct trfrSMSTraits {
	short	pr;
	short	prMethod; 
	short	memSize; 
	short goalType;
	float	dp;
	float	gb;
	float alphaDB;
	int betaDB; 
	float	stepMort;
	bool straightenPath;
};

// structures for settlement parameters

struct settleType {
	bool stgDep;
	bool sexDep; 
	bool indVar;
	short settTrait[2];
};
struct settleRules {
	bool densDep;
	bool wait; 
	bool goToNeighbourLocn; 
	bool findMate;
};
struct settleSteps {
	int minSteps;
	int maxSteps;
	int maxStepsYr;
};
struct settleTraits {
	float s0;
	float alpha; 
	float beta;

	settleTraits() : s0(0.0), alpha(0.0), beta(0.0) {};

	settleTraits(const settleTraits& e) : s0(e.s0), alpha(e.alpha), beta(e.beta) {};

	void divideTraitsBy(int i) {
		s0 /= i;
		alpha /= i;
		beta /= i;
	}
};

// Structures for interactions

// Initiated interactions,
// e.g. predation/parasitism/pollination 
// from P.O.V. of predator/parasite/pollinator
// Effector species "owns" all parameters of the functional response
struct initdInteraction {
	Species* recipientSpecies; // host, prey...
	int stage; // which recipient stage this applies to
	double beta; // conversion rate
	double handlingTime; // how many such interactions resolved per generation?
	double interactionRate; // e.g. attack rate
	double hullCoeff; // shape of the functional response
	double interfIntercept; 
	double interfExponent;
	double relPreference; // weight for choosing this prey over others
};

// Received interactions,
// e.g. predation/parasitism/pollination 
// from P.O.V. of prey/host/flowering plant
struct recdInteraction {
	Species* initiatorSpecies; // predator, pollinator, parasite...
	int stage; // which initiator stage this applies to
	double delta; // effect of one unit interaction on the process
};

// Resource-mediated interactions,
// e.g. scramble competition and mutualism
struct resInteraction {
	Species* partnerSpecies; // competitor or mutualist
	int stage; // which partner species stage this applies to
	double alpha; // effect of one individual of the partner species
};

//---------------------------------------------------------------------------
typedef short species_id;

class Species {

public:
	Species(const species_id& sp);
	~Species();
	species_id getID();

	// Demographic parameter functions

	// Create habitat carrying capacity table
	void createHabK(short nbHab);
	void setHabK(
		short habIx,	// may differ from habitat no. supplied by user
		float habK		// carrying capacity (inds/cell)
	);
	float getHabK(short habIx);
	float getMaxK(); // return highest carrying capacity over all habitats
	void deleteHabK();

	void setNbStages(const int nbStg) { nStages = nbStg; }
	void setStage(const stageParams stgPar);
	stageParams getStageParams();
	void setDemogr(const demogrParams dem);
	demogrParams getDemogrParams();
	short getRepType();
	bool stageStructured();
	void setDensDep(float devCoeff, float survCoeff);
	densDepParams getDensDep(); // Get development and survival coefficients

	void setFec(short stg, short sex, float fec); // stg must be > 0
	float getFec(short stg, short sex);
	void setDev(short stg, short sex, float devProb);
	float getDev(short stg, short sex);
	void setSurv(short stg, short sex, float survProb);
	float getSurv(short stg, short sex);
	float getMaxFec(); // Get highest fecundity of any stage
	void setMinAge(short stg, short sex, int minAge); // must be 0 for stages 0,1
	short getMinAge(short stg, short sex);
	void createDDwtFec(short dim);	// no. of stages * no. of sexes
	void setDDwtFec(short row, short col, float weight);
	float getDDwtFec(short row, short col);
	void deleteDDwtFec(); // Delete fecundity weights matrix
	void createDDwtDev(short dim);	// no. of stages * no. of sexes
	void setDDwtDev(short row, short col, float weight);
	float getDDwtDev(short row, short col);
	void deleteDDwtDev();
	void createDDwtSurv(short dim);	// no. of stages * no. of sexes
	void setDDwtSurv(short row, short col, float weight);
	float getDDwtSurv(short row, short col);
	void deleteDDwtSurv(); // Delete survival weights matrix

	// Environmental stochasticity
	void setMinMax(float min, float max);
	float getMinMax(short option); //0 = return minimum, otherwise = return maximum
	float getLocalExtProb() const { return localExtProb; }
	void setLocalExtProb(float p) { localExtProb = p; }

	// Patch sampling for genetics output
	string getSamplingOption() const { return patchSamplingOption; }
	bool doesOutputGeneValues() const { return output.outputGenes; }
	bool doesOutputWeirCockerham() const { return output.outputWeirCockerham; }
	bool doesOutputWeirHill() const { return output.outputWeirHill; }
	bool isGeneticOutputYear(int yr) const {
		if (output.outputGenes || output.outputWeirCockerham || output.outputWeirHill)
			return yr >= output.outputStartGenetics
			&& yr % output.outputGeneticInterval == 0;
		else return false;
	}
	std::set<int>& getSamplePatches() { return samplePatchList; }
	string getNIndsToSample() { return nIndsToSample; }
	std::set<int>& getStagesToSample() { return stagesToSampleFrom; }
	int getNbPatchesToSample() { return nPatchesToSample; }

	// Environmmental gradient
	void setEnvGrad(const envGradParams& envGrad) { grad = envGrad; }
	envGradParams getEnvGradient() const { return grad; }
	bool usesGradient() const { return grad.usesGradient; }
	bool isGradientShifting(int year) {
		return grad.doesShift && year > grad.shiftBegin && year < grad.shiftStop;
	}
	void incrementGradOptY() { grad.optY++; }
	void resetOptY() { grad.optY = grad.optY0; }

	// Genetic functions
	void resetGeneticParameters();
	bool areMutationsOn();
	bool isDiploid() const;
	int incrNbGenLoadTraits();
	int getNbGenLoadTraits() const;

	// Emigration parameter functions
	void setEmigRules(const emigRules emig);
	emigRules getEmigRules(); // Get emigration rules
	void setSpEmigTraits(short stg, short sex, const emigTraits emig);
	emigTraits getSpEmigTraits(short stg, short sex);
	float getSpEmigD0(short stg, short sex);

	// Transfer parameter functions
	void setTrfrRules(const transferRules trfr);
	transferRules getTransferRules();
	void setFullKernel(bool	usesFullKernel);
	bool useFullKernel();
	void setSpKernTraits(short stg, short sex, const trfrKernelParams trfr,	const int landResol);
	trfrKernelParams getSpKernTraits(short stg, short sex);
	void setMortParams(const trfrMortParams	trfrMort);
	trfrMortParams getMortParams();
	void setSpMovtTraits(const trfrMovtParams trfrMovement);
	trfrMovtParams getSpMovtTraits();
	trfrCRWTraits getSpCRWTraits();	
	trfrSMSTraits getSpSMSTraits();
	short getMovtHabDim(); 	// dimension of habitat-dependent step mortality and costs matrices
	void createHabCostMort(short nbHab);
	void setHabCost(short habIx, int cost);
	void setHabMort(short habIx, double habMortProb);
	int getHabCost(short habIx);
	double getHabMort(short habIx);
	void deleteHabCostMort();

	// Settlement parameter functions
	void setSettle(const settleType sett);
	settleType getSettle(); // Get settlement type
	void setSettRules(short stg, short sex, const settleRules sett);
	settleRules getSettRules(short stg, short sex);
	void setSteps(short stg, short sex, const settleSteps stepLimits);
	settleSteps getSteps(short stg, short sex);
	// Set settlement density dependence traits
	void setSpSettTraits(const short stg, const short sex, const settleTraits);
	// Get settlement density dependence traits
	settleTraits getSpSettTraits(short stg, short sex);

	// Genome management functions
	void addTrait(TraitType traitType, const SpeciesTrait& trait);
	void clearTraitTable();
	SpeciesTrait* getSpTrait(TraitType trait) const;
	std::set<TraitType> getTraitTypes();
	int getNTraits() const;
	int getNPositionsForTrait(const TraitType trait) const;
	int getGenomeSize() const;
	float getRecombinationRate() const;
	std::set<int> getChromosomeEnds() const;
	void setGeneticParameters(
		const std::set<int>& chromosomeEnds,
		const int genomeSize, 
		const float recombinationRate,
		string patchSamplingOption,
		const std::set<int>& samplePatchList, 
		const string nIndsToSample, 
		const std::set<int>& stagesToSampleFrom, 
		int nPatchesToSampleFrom
	);
	void setSamplePatchList(const std::set<int>& samplePatchList);

	// Output control functions
	bool doesOutputInds() const { return output.outInds; }
	bool isIndOutputYear(int yr) const {
		if (output.outInds) return yr >= output.outStartInd
			&& yr % output.outIntInd == 0;
		else return false;
	}
	bool doesOutputConnect() const { return output.outConnect; }
	bool isConnectOutputYear(int yr) const {
		if (output.outConnect) return yr >= output.outStartConn
			&& yr % output.outIntConn == 0;
		else return false;
	}
	bool savesVisits() const { return output.saveVisits;  }
	bool doesOutputOccup() const { return output.outOccup; }
	bool isOccOutputYear(int yr) const {
		if (output.outOccup) return yr % output.outIntOcc == 0;
		else return false;
	}
	int getOutOccInt() const { return output.outIntOcc; }
	bool doesOutputPop() const { return output.outPop; }
	bool isPopOutputYear(int yr) const { 
		if (output.outPop) return yr >= output.outStartPop
			&& yr % output.outIntPop == 0;
		else return false;
	}
	bool doesOutputRange() const { return output.outRange; }
	bool isRangeOutputYear(int yr) const {
		if (output.outRange) return yr % output.outIntRange == 0;
		else return false;
	}
	bool doesOutputTraitCell() const { return output.outTraitsCells; }
	bool isTraitCellOutYear(int yr) const {
		if (output.outTraitsCells) return yr >= output.outStartTraitCell
			&& yr % output.outIntTraitCell == 0;
		else return false;
	}
	int getOutTrCellInt() const { return output.outIntTraitCell; }
	bool doesOutputTraitRows() const { return output.outTraitsRows; }
	bool isTraitRowOutYear(int yr) const {
		if (output.outTraitsRows) return yr >= output.outStartTraitRow
			&& yr % output.outIntTraitRow == 0;
		else return false;
	}
	outputParams getOutputParams() const { return output; }
	void setOutputParams(const outputParams& out) { output = out; }
#if RS_RCPP
	bool doesOutputPaths() const { return output.outPaths; }
#endif

	void setInitParams(initParams initPars) { init = initPars; }
	initParams getInitParams() const { return init; }
	initInd getInitInd(int ix);
	void addInitInd(initInd iind);
	void resetInitInds();
	int getNbInitInds();
	void setProp(short stg, float initProp);
	float getProp(short stg);
	bool isRestrictYear(int yr) const {
		return init.restrictRange
			&& yr > init.initFrzYr
			&& yr < init.finalFrzYr;
	}
	void applyRangeRestriction();
	void resetRangeRestrictions(int dimX, int dimY);
	void freezeYrange(int minYlimit, int maxYlimit, int maxXlimit);
	bool isWithinLimits(const int& x, const int& y);

	// Interaction functions
	void addInteractingSpecies(const species_id& sp);
	set<species_id> getInteractingSpecies() const;

private:
	species_id ID;

	// Demographic parameters
	short repType;		// 0 = asexual, 1 = simple two sex, 2 = complex two sex
	short nStages;      // no. of stages (incl. juveniles) in structured population
	float propMales;    // proportion of males at birth in sexual model
	float harem;        // max harem size in complex sexual model
	float bc;			// competition coefficient for non-structured population
	float lambda;       // max intrinsic growth rate for non-structured population
	float probRep; 		// probability of reproducing in subsequent seasons
	short repSeasons;	// no. of reproductive seasons per year
	short repInterval;	// no. of reproductive seasons between subsequent reproductions
	short maxAge;       // max age in structured population
	short survival;		// survival timing: 0 = at reprodn, 1 = between reprodns, 2 = anually
	bool stageStruct;
	bool fecDens;
	bool fecStageDens;
	bool devDens;
	bool devStageDens;
	bool survDens;
	bool survStageDens;
	bool disperseOnLoss;	// individuals disperse on complete loss of patch
							// (otherwise they die)
	short habDimK;		// dimension of carrying capacities matrix
	float* habK;		// habitat-specific carrying capacities (inds/cell)
	float devCoeff; 	// density-dependent development coefficient
	float survCoeff; 	// density-dependent survival coefficient
	float** ddwtFec;    // density-dependent weights matrix for fecundity
	float** ddwtDev;    // density-dependent weights matrix for development
	float** ddwtSurv;   // density-dependent weights matrix for survival
	// NB for the following arrays, sex 0 is females, sex 1 is males
	float fec[gMaxNbStages][gMaxNbSexes];		// fecundities
	float dev[gMaxNbStages][gMaxNbSexes];		// development probabilities
	float surv[gMaxNbStages][gMaxNbSexes];		// survival probabilities
	short minAge[gMaxNbStages][gMaxNbSexes];	// minimum age to enter stage
	short ddwtFecDim;	// dimension of density-dependent weights matrix for fecundity
	short ddwtDevDim;	// dimension of density-dependent weights matrix for fecundity
	short ddwtSurvDim;	// dimension of density-dependent weights matrix for fecundity
	float minRK; 		// minimum growth rate OR carrying capacity
	float maxRK; 		// maximum (under environmental stochasticity)
	float localExtProb; // probability of any population to go extinct on a given year
	envGradParams grad;

	// Genome parameters
	/** The traits table. */
	std::map<TraitType, std::unique_ptr<SpeciesTrait>> spTraitTable;
	std::set<int> chromosomeEnds;
	int genomeSize;
	bool diploid;
	bool mutationsOn;
	int nbGeneticFitnessTraits;
	float recombinationRate;
	string patchSamplingOption;
	std::set<int> samplePatchList;
	int nPatchesToSample; // for cell-based landscape
	std::set<int> stagesToSampleFrom;
	string nIndsToSample; // could be integer or 'all', all means in in selected patches not necessarily all in population

	// Emigration parameters
	bool densDepEmig;	// density-dependent emigration
	bool stgDepEmig;	// stage-dependent emigration
	bool sexDepEmig;	// sex-dependent emigration
	bool indVarEmig;	// individual variation in emigration
	short emigStage;	// stage which emigrates (used for stage-strucutred population
						// having individual variability in emigration probability)
	// sex 0 is females, sex 1 is males
	float	d0[gMaxNbStages][gMaxNbSexes];			 // maximum emigration probability
	float	alphaEmig[gMaxNbStages][gMaxNbSexes];	 // slope of density-dependent reaction norm
	float	betaEmig[gMaxNbStages][gMaxNbSexes];	 // inflection point of reaction norm (in terms of N/K)
	
	// Transfer parameters
	bool usesMovtProcess;
	bool stgDepTrfr;
	bool sexDepTrfr;
	bool distMort;
	bool indVarTrfr;
	bool twinKern;
	bool habMort;		// habitat-dependent mortality
	float	meanDist1[gMaxNbStages][gMaxNbSexes];	// mean of 1st dispersal kernel (m)
	float	meanDist2[gMaxNbStages][gMaxNbSexes];	// mean of 2nd dispersal kernel (m)
	float	probKern1[gMaxNbStages][gMaxNbSexes];	// probability of dispersing with the 1st kernel
	// As evolving traits are are not stage-dependent, no. of rows can be 1

	float fixedMort;		// constant mortality probability
	float mortAlpha;		// slope for mortality distance dependence function
	float mortBeta;			// inflection point for mortality distance dependence function
	short moveType; 		// 1 = SMS, 2 = CRW
	short pr;				// SMS perceptual range (cells)
	short prMethod;			// SMS perceptual range evaluation method:
							// 1 = arith. mean, 2 = harmonic mean, 3 = inverse weighted arith. mean
	short memSize;			// SMS memory size (1-14 steps)
	short goalType;			// SMS goal bias type: 0 = none, 1 = towards goal, 2 = dispersal bias
	float dp;				// SMS directional persistence
	float gb;				// SMS goal bias
	float alphaDB; 			// SMS dispersal bias decay rate
	int betaDB; 			// SMS dispersal bias decay inflection point (no. of steps)
	float stepMort;			// constant per-step mortality probability for movement models
	double* habStepMort;	// habitat-dependent per-step mortality probability
	float stepLength;		// CRW step length (m)
	float rho;				// CRW correlation coefficient
	short habDimTrfr;		// dimension of habitat-dependent step mortality and costs matrices
	int* habCost;			// habitat costs
	bool usesCosts;			// import cost map from file?
	bool straightenPath;	// straighten path after decision not to settle
	bool fullKernel;		// used to indicate special case when density-independent emigration
	// is 1.0, and kernel-based movement within the natal cell is used
	// to determine philopatry

	// Settlement parameters
	bool stgDepSett;
	bool sexDepSett;
	bool indVarSett;   									// individual variation in settlement
	bool densDepSett[gMaxNbStages][gMaxNbSexes];
	bool wait[gMaxNbStages][gMaxNbSexes];				// wait to continue moving next season (stage-structured model only)
	bool goToNeighbourLocn[gMaxNbStages][gMaxNbSexes];	// settle in neighbouring cell/patch if available (ditto)
	bool findMate[gMaxNbStages][gMaxNbSexes];
	int minSteps[gMaxNbStages][gMaxNbSexes];     		// minimum no. of steps
	int maxSteps[gMaxNbStages][gMaxNbSexes];			// maximum total no. of steps
	int maxStepsYr[gMaxNbStages][gMaxNbSexes]; 			// maximum no. of steps in any one dispersal period
	float s0[gMaxNbStages][gMaxNbSexes];				// maximum settlement probability
	float alphaS[gMaxNbStages][gMaxNbSexes];			// slope of the settlement reaction norm to density
	float betaS[gMaxNbStages][gMaxNbSexes];				// inflection point of the settlement reaction norm to density

	// Interaction parameters
	vector<initdInteraction> initiatedIntrcts;	// species targeting this, e.g. predators and pollinators
	vector<recdInteraction> recipientIntrcts;	// species this one targets, e.g. preys or hosts
	vector<resInteraction> resDepIntrcts;	// species that affect carrying capacity (competitors/mutualists)

	// Initialisation parameters
	initParams init;
	vector<initInd> initInds;	// individuals to be initialised
	vector<float> initProps; // proportion of each stage among initial individuals
	int minX, minY, maxX, maxY; // limits the species is allowed to live in,
								// changed by the freeze/range restriction feature
	// Output controls
	outputParams output;
};

// Map to record and track all species
typedef map<species_id, Species*> speciesMap_t;

//---------------------------------------------------------------------------

#ifndef NDEBUG
// For testing purposes only
demogrParams createDefaultHaploidDemogrParams();
demogrParams createDefaultDiploidDemogrParams();
#endif // NDEBUG

//---------------------------------------------------------------------------
#endif
