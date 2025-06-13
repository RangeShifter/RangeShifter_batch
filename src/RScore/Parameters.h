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

 RangeShifter v2.0 Parameters

 Implements the following classes:

 paramSim   - Simulation parameters
 paramStoch - Environmental stochasticity parameters

 Also declares some structures and functions used throughout the program.

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 25 June 2021 by Steve Palmer

 ------------------------------------------------------------------------------*/

#ifndef ParametersH
#define ParametersH

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <algorithm>
using namespace std;

#include "RSrandom.h"

constexpr int gNoDataCost = 100000; // cost to use in place of nodata value for SMS;
constexpr int gAbsorbingNoDataCost = 100; // cost to use in place of nodata value for SMS;
// when boundaries are absorbing
constexpr int gMaxNbStages = 10;		// maximum number of stages permitted
constexpr int gMaxNbSexes = 2;			// maximum number of sexes permitted

#if RS_RCPP
#ifndef R_EXT_CONSTANTS_H_  // the R headers define PI as a macro, so that the 'else' line results in an error
#define M_2PI 6.283185307179586
const double PI = 3.141592653589793238462643383279502884197169399375;
#endif
#else
#define M_2PI 6.283185307179586
const double PI = 3.141592654;
#endif

const double SQRT2 = std::sqrt(double(2.0)); // more efficient than calculating every time

//---------------------------------------------------------------------------

// Common declarations
struct locn { int x; int y; };

//--------------------------------------------------------------------------

/** Trait types **/

enum TraitType {
	NEUTRAL, 
	GENETIC_LOAD, GENETIC_LOAD1, GENETIC_LOAD2, GENETIC_LOAD3, GENETIC_LOAD4, GENETIC_LOAD5,

	E_D0, E_ALPHA, E_BETA,
	S_S0, S_ALPHA, S_BETA,

	E_D0_F, E_ALPHA_F, E_BETA_F,
	S_S0_F, S_ALPHA_F, S_BETA_F,

	E_D0_M, E_ALPHA_M, E_BETA_M,
	S_S0_M, S_ALPHA_M, S_BETA_M,

	CRW_STEPLENGTH, CRW_STEPCORRELATION,

	KERNEL_MEANDIST_1, KERNEL_MEANDIST_2, KERNEL_PROBABILITY,
	KERNEL_MEANDIST_1_F, KERNEL_MEANDIST_2_F, KERNEL_PROBABILITY_F,
	KERNEL_MEANDIST_1_M, KERNEL_MEANDIST_2_M, KERNEL_PROBABILITY_M,

	SMS_DP, SMS_GB, SMS_ALPHADB, SMS_BETADB,

	INVALID_TRAIT // error
};

enum GenParamType { MEAN, SD, MIN, MAX, SHAPE, SCALE, INVALID };
enum DistributionType { UNIFORM, NORMAL, GAMMA, NEGEXP, SCALED, KAM, SSM, NONE };
enum ExpressionType { AVERAGE, ADDITIVE, NOTEXPR, MULTIPLICATIVE };

string to_string(const TraitType& tr);
string to_string(const GenParamType& param);
string to_string(const DistributionType& dist);
string to_string(const ExpressionType& expr);

/** Param's types **/
typedef enum { KERNEL, SMS, CRW} movement_t;

//sex types
typedef enum {
	FEM = 0, MAL = 1,
	NA, // not applicable. e.g. for NEUTRAL or genetic load trait
	INVALID_SEX // error
} sex_t;

//---------------------------------------------------------------------------

// Environmental gradient parameters

struct envGradParams {
	bool usesGradient; 
	bool doesShift;
	int gradType;	// 0 = none
					// 1 = carrying capacity
					// 2 = growth rate
					// 3 = local extinction probability
	float gradIncr; // gradient steepness
	float optY;		// current optimum row
	float optY0;	// optimum row before shift
	float factor;	// local scaling factor
	float extProbOpt;
	float shiftRate; // rows per year
	int shiftBegin; 
	int shiftStop;
};

//----------------------------------------

// Environmental stochasticity parameters

struct envStochParams {
	bool usesStoch; 
	bool stochIsLocal;
	bool inK; 
	float ac; 
	float std;
};

class paramStoch {

public:
	paramStoch();
	~paramStoch();
	void setStoch(envStochParams);
	bool envStoch();
	envStochParams getStoch();

private:
	bool usesStoch;		// stochasticity applied
	bool stochIsLocal;	// applied locally (if not, application is global)
	bool stochInK;		// in carrying capacity (if not, in growth rate)
	float ac;			// temporal autocorrelation coefficient		
	float std;			// amplitude of fluctuations: sampled from N(0,std)
};

//---------------------------------------------------------------------------

// Initialisation (seeding) parameters

struct initParams {
	short seedType;		// initialisation type: 
	//	0 = free, 
	//	1 = from species distn,
	//	2 = initial individuals, 
	//	3 = from file
	short freeType;		// free initialisation type:
	//	0 = random (given no.)
	//	1 = all suitable cells/patches
	//	2 = manually selected cells/patches
	short spDistType;	// species distribution initialisation type:
	//	0 = all suitable cells/patches,
	//	1 = some randomly chosen suitable cells/patches,
	//	2 = all cells/patches within selected sp. dist. cells
	short initDens;		// initialisation density:
	//	0 = at carrying capacity
	//	1 = at half carrying capacity
	//	2 = specified no. per cell or density
	short initAge;		// initial age distribution within each stage:
	//	0 = lowest possible age
	//	1 = randomised
	//	2 = quasi-equilibrium
	int initFrzYr;		 	// year until which initial range is frozen
	bool restrictRange;		// restrict range to northern front
	int restrictRows;		// no. of rows to retain behind front
	int restrictFreq;		// frequency of range restriction
	int finalFrzYr;		 	// year after which range is frozen
	int indsCell;			 // initial individuals / cell (cell-based model)
	float indsHa;			 // initial density (patch-based model)
	int minSeedX;			 // min. and max. of area to initialise (cell numbers),
	int maxSeedX;			 // only applied if seedType is 0
	int minSeedY;
	int maxSeedY;
	int nSeedPatches;	 	// no. of cells/patches to initialise
	int nSpDistPatches;		// no. of species distribution cells to initialise
	string indsFile;		// no. of species distribution cells to initialise
	
};

struct initInd {
	int year, patchID, x, y; 
	short sex, age, stage;
	short speciesID;
};

//---------------------------------------------------------------------------

// Simulation parameters

struct simParams {
	int batchNum;
	int simulation; 
	int reps; 
	int years;
	bool usesStageStruct;
	bool absorbing;
	bool fixReplicateSeed;
	bool batchMode;
};

class paramSim {

public:
	paramSim(const string& pathToProjDir = "");
	~paramSim();
	void setSim(simParams sim);
	simParams getSim();
	int getSimNum();
	string getDir(int option);
	void setBatchNum(const int& batchNb) {
		batchNum = batchNb;
		batchMode = true;
	}
#if RS_RCPP
	bool getReturnPopRaster();
	bool getCreatePopFile();
#endif

private:
	int batchNum;				// batch number
	int simulation;				// simulation no.
	int reps;					// no. of replicates
	int years;					// no. of years
	bool batchMode;		
	bool usesStageStruct;       // is stage-structure enabled?
	bool absorbing; 			// landscape boundary and no-data regions are absorbing boundaries
	string dir;					// full name of working directory
	bool fixReplicateSeed;
};

struct outputParams {
	int outStartPop;			// output start year for population file
	int outStartInd;			// output start year for individuals file
	int outStartTraitCell;		// output start year for traits by cell file
	int outStartTraitRow;		// output start year for traits by row file
	int outStartConn;			// output start year for connectivity matrix
	int outIntRange;			// output interval for range file
	int outIntOcc;				// output interval for occupancy file
	int outIntPop;				// output interval for population file
	int outIntInd;				// output interval for individuals file
	int outIntTraitCell;		// output interval for traits by cell file
	int outIntTraitRow;			// output interval for traits by row file
	int outIntConn;				// output interval for connectivity matrix
	int traitInt;				// output interval for evolving traits maps
	bool outRange;				// produce output range file?
	bool outOccup;				// produce output occupancy file?
	bool outPop;				// produce output population file?
	bool outInds;				// produce output individuals file?
	bool outTraitsCells;		// produce output summary traits by cell file?
	bool outTraitsRows;			// produce output summary traits by row (y) file?
	bool outConnect;			// produce output connectivity file?
	bool saveVisits;			// save dispersal visits heat maps?
	bool outputGenes;
	bool outputWeirCockerham;
	bool outputWeirHill;
	int outputStartGenetics;
	int outputGeneticInterval;
#if RS_RCPP
	int outStartPaths;
	int outIntPaths;
	bool outPaths;
	bool ReturnPopRaster;
	bool CreatePopFile;
#endif
};

extern RSrandom* pRandom;

//---------------------------------------------------------------------------
#endif
