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

 RangeShifter v2.0 BatchMode

 Functions for running in BATCH MODE

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 28 July 2021 by Greta Bocedi

 ------------------------------------------------------------------------------*/

#ifndef BatchModeH
#define BatchModeH

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <regex>
using namespace std;

#include "./RScore/Parameters.h"
#include "./RScore/Landscape.h"
#include "./RScore/Species.h"
#include "./RScore/Model.h"
#include "./RScore/SpeciesTrait.h"
#include "./RScore/NeutralTrait.h"

// Global variables to check parameter
// consistency across input files
struct spInputOptions {

	int reproType;
	int nbStages;

	// Track trait-relevant options to check for coherency across input files, 
	// e.g. if emig file says emigration is indvar, trait file should have d0 entry
	bool anyNeutral = false;
	bool isEmigIndVar = false;
	bool isEmigDensDep = false;
	bool isEmigSexDep = false;
	bool isSettIndVar = false;
	bool isSettSexDep = false;
	bool isKernTransfIndVar = false;
	bool isKernTransfSexDep = false;
	bool usesTwoKernels = false;
	bool isSMSTransfIndVar = false;
	bool usesSMSGoalBias = false;
	bool isCRWTransfIndVar = false;
	int nbTraitFileRows; // how many lines to expect in traitsFile?
};

bool traitExists(const TraitType& tr, const vector<TraitType>& existingTraits);
TraitType addSexDepToTrait(const TraitType& t, const sex_t& sex);
int checkTraitSetCoherency(const vector <TraitType>& allReadTraits, const int& simNb, const species_id& sp);

constexpr int gEmptyVal = -9;

struct simCheck {
	bool isNewSim;
	int simNb, simLines, reqdSimLines, errors;
};

bool checkInputFiles(string pathToControlFile, string inputDir, string outputDir);
bool CheckSimFile();
bool CheckParameterFile();
bool CheckLandFile(int landType, string inputDir);
bool CheckSpLandFile(string inputDir, bool isInitial);
int CheckGeneticsFile(string inputDir);
int CheckDynamicFile(string inputDir);
int CheckStageFile(string inputDir);
bool CheckTransitionFile(short, short nbSexesDemogr);
bool CheckWeightsFile(string fileType, int nbStages, int nbSexes);
int CheckEmigFile();
int CheckTransferFile(string inputDir);
int CheckSettleFile();
int CheckInitFile(string inputDir);
int CheckInitIndsFile(int simNb, species_id sp);
simCheck CheckStageSex(string, int, int, species_id sp, simCheck, int, int, int, int, int, bool, bool);
int CheckGeneticsFile(string inputDir);
int CheckTraitsFile(string inputDir);

void BatchError(
	string,	// file name
	int,		// line number
	int,		// option
	string	// fieldname
);
/* Write error message to batch log file. Options are as follows:
0 - general message only, no reference to field name
1 - fieldname must be 0 or 1
2 - fieldname must be 0, 1 or 2
3 - fieldname must be 0, 1, 2 or 3
4 - fieldname must be from 0 to 4
5 - fieldname must be from 0 to 5
6 - fieldname must be from 0 to 6
7 - fieldname must be from 0 to 7
10 - fieldname must be >0
11 - fieldname must be >=1
12 - fieldname must be >=2
13 - fieldname must be >=3
18 - fieldname must be >1
19 - fieldname must be >=0
20 - fieldname must be between 0 and 1
21 - fieldname must be >1
33 - fieldname must be 1, 2 or 3
44 - fieldname must be from 1 to 4
55 - fieldname must be from 1 to 5
66 - fieldname must be from 1 to 6
100 - fieldname must be between 0 and 100
111 - simulation must match first simulation in ParameterFile
222 - simulations must be sequential integers
223 - seasons must be sequential integers
333 - columns must match no. of habitats
444 - columns must be one fewer than no. of stages
555 - columns must match no. of stages
666 - fieldname must be a unique positive integer
*/
void BatchError(
	string,	// file name
	int,		// line number
	int,		// option
	string,	// fieldname
	string	// fieldname2
);
/* Write error message to batch log file. Options are as follows:
1 - fieldname must be greater than fieldname2
2 - fieldname must be greater than or equal to fieldname2
3 - fieldname must be less than or equal to fieldname2
4 - fieldname must be less than fieldname2
*/

bool isValidFractalDim(int x);

void printControlFormatError();
void FormatError(string, int);
void OpenError(string, string);
void EOFerror(string);
void FileOK(string, int, int);
void FileHeadersOK(string);
void SimulnCountError(string);

void RunBatch();
int ReadParameters(const Landscape*, speciesMap_t& simSpecies);
int ReadLandFile(Landscape*);
void ReadSpLandFile(
	ifstream& ifsSpLand,
	map<species_id, string>& pathsToPatchMaps,
	map<species_id, string>& pathsToCostMaps,
	map<species_id, string>& pathsToSpDistMaps,
	map<species_id, bool>& whichUseSpDist
);
int ReadDynLandFile(Landscape*);
int ReadStageStructure(speciesMap_t& simSpecies);
int ReadTransitionMatrix(
	Species* pSpecies,
	short nstages, 
	short nsexesDem, 
	short hab, 
	short season
);
int ReadStageWeights(Species* pSpecies, int option);
int ReadEmigration(speciesMap_t& simSpecies);
int ReadTransferFile(speciesMap_t& simSpecies, landParams paramsLand, int transferType);
int ReadTransferKernels(speciesMap_t& simSpecies, landParams paramsLand);
void ReadTransferSMS(speciesMap_t& simSpecies, const landParams&);
int ReadTransferCRW(speciesMap_t& simSpecies, const landParams&);
int ReadSettlement(speciesMap_t& simSpecies);
int ReadInitialisation(const landParams& paramsLand, speciesMap_t& simSpecies);
int ReadInitIndsFile(
	Species* pSpecies, 
	int option, 
	const landParams& paramsLand,
	string indsfile
);
int ReadGeneticsFile(speciesMap_t& simSpecies, ifstream& ifs);
int ReadTraitsFile(speciesMap_t& simSpecies, ifstream& ifs, const int& whichSim);

// Helper functions to ReadGenetics and ReadTraits
void setUpSpeciesTrait(Species* pSpecies, vector<string>);
DistributionType stringToDistributionType(const std::string& str);
ExpressionType stringToExpressionType(const std::string& str);
map<GenParamType, float> stringToParameterMap(string parameters);
set<int> selectRandomLociPositions(int noLoci, const int& genomeSize);
set<int> stringToLoci(string pos, string nLoci, const int& genomeSize);
TraitType stringToTraitType(const std::string& str);
const sex_t stringToSex(const std::string& str);
set<int> stringToPatches(const string&);
set<int> stringToStages(const string&, const int&);
set<int> stringToChromosomeEnds(string, const int&);
GenParamType strToGenParamType(const string& str);

// external pointers to parameter sets
extern paramStoch* paramsStoch;
extern paramSim* paramsSim;
extern int RS_random_seed;

//---------------------------------------------------------------------------
#endif
