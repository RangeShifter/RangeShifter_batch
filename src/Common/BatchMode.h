/*------------------------------------------------------------------------------

RangeShifter v2.0 BatchMode

Functions for running in BATCH MODE

For full details of RangeShifter, please see:
Bocedi G., Palmer S.C.F., Pe�er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species� responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 22 August 2018 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef BatchModeH
#define BatchModeH

#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
using namespace std;

#include "Parameters.h"
#include "Landscape.h"
#include "Species.h"
#include "Model.h"
#if RS_ABC
#include "ABC.h"
#endif

struct batchfiles {
	bool ok;
	int batchNum;
	int nSimuls;
	int nLandscapes;
	int patchmodel,resolution,landtype,maxNhab,speciesdist,distresolution;
	int reproductn,repseasons,stagestruct,stages,transfer;
	int sexesDem;		// no. of explicit sexes for demographic model
	int sexesDisp;	// no. of explicit sexes for dispersal model
	string parameterFile;
	string landFile;
	string stageStructFile;
	string emigrationFile;
	string transferFile;
	string settleFile;
	string geneticsFile;
	string initFile;
#if VIRTUALECOLOGIST
	string virtEcolFile;
#endif
#if RS_ABC
	string abcParamsFile,abcObsFile;
#endif
};

struct simCheck {
	bool newsimul;
	int simul,simlines,reqdsimlines,errors;
};

batchfiles ParseControlFile(string,string,string);
#if BUTTERFLYDISP
int ParseParameterFile(string);
#else
int ParseParameterFile(void);
#endif
int ParseLandFile(int,string);
int ParseDynamicFile(string);
int ParseStageFile(string);
int ParseTransitionFile(void);
int ParseWeightsFile(string);
int ParseEmigFile(void);
int ParseTransferFile(string);
int ParseSettleFile(void);
int ParseGeneticsFile(string);
int ParseArchFile(void);
int ParseInitFile(string);
int ParseInitIndsFile(void);
#if VIRTUALECOLOGIST
int ParseVirtEcolFile(string);
int ParseSampleFile(void);
int ParsePatchFile(void);
#endif // VIRTUALECOLOGIST
#if RS_ABC
int ParseABCParamsFile(void);
int ParseABCObsFile(void);
#endif // RS_ABC
#if EVOLSMS
int ParseMortFile(void);
int ReadMortalities(string);
#endif // EVOLSMS
int CheckCostRaster(string,string);
simCheck CheckStageSex(string,int,int,simCheck,int,int,int,int,int,bool);

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
333 - columns must match no. of habitats
444 - columns must match no. of stages
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

int power2check(int x);

void CtrlFormatError(void);
void ArchFormatError(void);
#if VIRTUALECOLOGIST
void SampleFormatError(void);
#endif // VIRTUALECOLOGIST
void FormatError(string,int);
void OpenError(string,string);
void EOFerror(string);
void FileOK(string,int,int);
void FileHeadersOK(string);
void SimulnCountError(string);

void RunBatch(int,int);
int ReadParameters(int,Landscape*);
int ReadLandFile(int);
int ReadLandFile(int,Landscape*);
int ReadDynLandFile(Landscape*);
int ReadStageStructure(int);
int ReadTransitionMatrix(int);
int ReadStageWeights(int);
int ReadEmigration(int);
int ReadTransfer(int,Landscape*);
int ReadSettlement(int);
int ReadGenetics(int);
int ReadArchFile(string);
int ReadInitialisation(int,Landscape*);
int ReadInitIndsFile(int,Landscape*,string);
#if VIRTUALECOLOGIST
int ReadVirtEcol(int);
int ReadSampleFile(string);
int ReadPatchFile(string);
#endif // VIRTUALECOLOGIST

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

// external pointers to parameter sets
extern paramGrad *paramsGrad;
extern paramStoch *paramsStoch;
extern paramInit *paramsInit;
extern paramSim *paramsSim;

extern Species *pSpecies;
extern string costmapname;	// see FormMove.cpp (VCL) OR Main.cpp (batch)
extern string genfilename;	// see FormGenetics.cpp (VCL) OR Main.cpp (batch)
#if VIRTUALECOLOGIST
extern string locfilename;		// see FormVirtEcol.cpp (VCL) OR Main.cpp (batch)
extern string patchfilename;	// see [NOT YET CODED FOR GUI] (VCL) OR Main.cpp (batch)
#endif // VIRTUALECOLOGIST
#if EVOLSMS
extern string mortfilename;	// see [NOT YET CODED FOR GUI] (VCL) OR Main.cpp (batch)
#endif // EVOLSMS
#if GROUPDISP || RS_ABC
#if RSRANDOM
extern int RS_random_seed;			// see RSrandom.cpp 
#endif
#endif // GROUPDISP

//---------------------------------------------------------------------------
#endif