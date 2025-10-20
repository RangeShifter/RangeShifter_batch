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

 RangeShifter v2.0 Landscape

 Implements the following classes:

 InitDist  - Initial species distribution

 Landscape - Landscape grid

 The Landscape is a rectangular array of Cells which are grouped together in
 Patches. As far as the user is aware, the Landscape is either patch-based or
 cell-based (having no Patches), but internally Patches are always present (they
 each comprise only one Cell for a cell-based model). The grain of the Landscape
 may be any positive integer, and is nominally in metres.

 The Landscape is either input from one or more text files in ArcGIS raster export
 format, or it is generated artificially as a random or fractal binary array (in
 which case, it must be cell-based). An input 'real' Landscape may hold within each
 Cell either discrete habitat classes, or percent cover of habitat classes, or a
 continuous quality index (1 to 100%).

 The Landscape may be dynamic, in which case the user specifies a set of changes
 to be applied at certain years during each simulation. The changes may be to
 habitat only, patches only (if a patch-based model) or habitats and patches.
 Although the changes must be supplied as entire habitat / patch files (which
 must match the original Landscape in terms of cell size and extent), internally
 they are recorded as lists of dynamic habitat and patch changes.

 The initial species distribution is a rectangular array if distribution cells
 (DistCell) covering the same spatial extent at the Landscape. Its grain may be
 either that of the Landscape or an integer multiple of it.

 The Landscape can record a list (in the vector initcells) of Cells or Patches
 to be intialised, which are specified by the user in FormSeeding. This option is
 available in the GUI version only.

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 28 July 2021 by Greta Bocedi
 ------------------------------------------------------------------------------*/

#ifndef LandscapeH
#define LandscapeH

#include <algorithm>
#include <fstream>
#include <vector>
#include <ranges>

using namespace std;

#include "Parameters.h"
#include "Patch.h"
#include "Cell.h"
#include "Species.h"
#include "FractalGenerator.h"
#if RS_RCPP
#include <locale>
#if !RSWIN64
#include <codecvt>
#endif
#include <Rcpp.h>
#endif

constexpr species_id gSingleSpeciesID = 0;

//---------------------------------------------------------------------------

// Initial species distribution

class InitDist {
public:
	InitDist();
	~InitDist();

	int readDistribution(string filename);
	void setDistribution(int nbInitialCells); // (0 for all cells)
	int cellCount();

	// Return the co-ordinates of a specified initial distribution
	// cell if it has been selected
	// otherwise return negative co-ordinates
	locn getSelectedCell(int cellIndex);
	int getResol() const { return resol; }

private:
	int resol;			// species distribution cell size (m)
	int maxX, maxY;		// dimensions
	double minEast;		// real world min co-ordinates
	double minNorth;	// read from raster file

	// list of cells in the initial distribution
	// cells MUST be loaded in the sequence ascending x within descending y
	std::vector <DistCell*> cells;
};

//---------------------------------------------------------------------------

struct landParams {
	bool usesPatches; 
	bool isArtificial;
	bool isDynamic;
	int landNum; 
	int resol;
	int nHab; 
	int nHabMax;
	int dimX, dimY, minX, minY, maxX, maxY;
	short rasterType;
};
struct landData {
	int resol; 
	int dimX, dimY, minX, minY, maxX, maxY;
};

bool isInLandBounds(const int& x, const int& y, const landData& land);

struct genLandParams {
	bool isFractal; bool isContinuous;
	float minPct, maxPct; float propSuit; float hurst; int maxCells;
};
struct landOrigin {
	double minEast; double minNorth;
};
struct rasterdata {
	bool ok;
	int errors, ncols, nrows, cellsize;
	double xllcorner, yllcorner;
#if RS_RCPP
	bool utf;
#endif
};
struct patchData {
	Patch* pPatch; 
	int patchNum, nCells; 
	int x, y;
};

struct cellChange {
	short originVal, currentVal, nextVal;
};

struct landChange {
	int chgnum = 0, chgyear = 999999;
	string habfile, spLandFile;
};
struct patchChange {
	int chgNb, x, y, oldPatch, newPatch;
};
struct costChange {
	int chgnum, x, y, oldcost, newcost;
};

class Landscape {
public:
	Landscape(const set<species_id>& speciesNames);
	~Landscape();
	void resetLand();

	// Generate patches, sample patches and set landscape limits
	void initialise(speciesMap_t& allSpecies, landParams land);
	void calcPatchOverlap();

	// Landscape parameters
	void setLandParams(landParams ppp, bool batchmode);
	landParams getLandParams() const;
	landData getLandData();
	void setGenLandParams(genLandParams ppp);
	genLandParams getGenLandParams() const;
	landOrigin getOrigin();

	// Habitat codes
	bool habitatsIndexed();
	void listHabCodes();
	void addHabCode(int hab);
	int findHabCode(int hab);
	int getHabCode(int xhab);
	void clearHabitats();

	// Patches and cells
	void setCellArray();
	void generatePatches(const speciesMap_t& allSpecies); // create patches for an artificial landscape
	void allocatePatches(const speciesMap_t& allSpecies); // create patches for a cell-based landscape
	Patch* addNewPatch(species_id id, int num);
	Patch* addNewPatch(species_id id, int seqnum, int num);
	void resetPatchLimits(species_id sp);
	void addNewCellToLand(int x, int y, float habQual);
	void addNewCellToLand(int x, int y, int habType);
	void addCellToLand(Cell* pCell);
	void addCellToPatch(species_id whichSpecies, Cell* pCell, Patch* pPatch);
	void addCellToPatch(species_id whichSpecies, Cell* pCell, Patch* pPatch, float habQual);
	void addCellToPatch(species_id whichSpecies, Cell* pCell, Patch* pPatch, int habType);
	void addNewCellToPatch(Patch* pPatch, int x, int y, int habType);
	void addNewCellToPatch(Patch* pPatch, int x, int y, float habQual);
	patchData getPatchData(species_id id, int patchIx);
	bool existsPatch(species_id whichSpecies, int patchIx);
	Patch* findPatch(species_id whichSpecies, int patchIx);
	set<int> getPatchNbs(species_id sp) const;
	void samplePatches(Species* pSpecies);
	int checkTotalCover();
	void resetPatchPopns();
	void updateCarryingCapacity(Species* pSpecies, int year, short landIx);
	Cell* findCell(int x, int y);
	int getPatchCount(species_id id) const;
	int allPatchCount() const;
	void updateHabitatIndices();

	// Environmental gradient
	void drawGradientDev();
	void updateEnvGradient(Species* pSpecies);

	// Environmental stochasticity
	void setGlobalStoch(int	nbYears);
	float getGlobalStoch(int year);
	void updateLocalStoch();

	// SMS costs
	void resetCosts();
	void resetEffCosts(species_id sp);

	// Dynamic changes
	void setDynamicLand(bool isDynamic);
	void addLandChange(landChange c);
	int numLandChanges();
	landChange getLandChange(short changeIx);
	void deleteLandChanges();
#if RS_RCPP && !R_CMD
	int readLandChange(
		int,		// change file number
		bool,		// change SMS costs?
		wifstream&, // habitat file stream
		wifstream&, // patch file stream
		wifstream&, // cost file stream
		int,		// habnodata
		int,		// pchnodata
		int			// costnodata
	);
#else
	int readLandChange(int fileNb, bool changeCosts);
#endif
	void createPatchChgMatrix();
	void resetPatchChanges();
	void recordPatchChanges(int landIx);
	int getNbPatchChanges(species_id sp);
	patchChange getPatchChange(species_id sp, int changeIx);
	void applyPatchChanges(species_id sp, const int& landChgNb, int& iPatchChg);
	void createCostsChgMatrix();
	void resetCostChanges();
	void recordCostChanges(int landIx);
	int getNbCostChanges(species_id sp);
	costChange getCostChange(species_id sp, int i);
	void applyCostChanges(species_id sp, const int& landChgNb, int& iCostChg);

	// Species distributions
	int newDistribution(species_id sp, string distFileName);
	void setDistribution(species_id sp, int nInit);
	// Specified cell matches one of the distn cells to be initialised?
	int distnCount();	// Return no. of initial distributions in the Landscape
	// Return no. of distribution cells in an initial distribution
	int distCellCount(species_id sp);
	// Get co-ordinates of a specified cell in a specified initial distn
	// Returns negative co-ordinates if the cell is not selected
	locn getSelectedDistnCell(species_id sp, int ix);
	bool usesSpDist(species_id sp) const;
	int getSpDistResol(species_id sp) const { return distns.at(sp).getResol(); }

	// Functions to handle connectivity matrix
	void createConnectMatrix(species_id sp);
	void resetConnectMatrix();
	void incrConnectMatrix(const species_id& speciesID, int originPatchNb, int settlePatchNb);
	void deleteConnectMatrix(const species_id& id);
	void outConnectHeaders(species_id sp);
	bool closeConnectOfs(species_id sp);
	void outConnect(species_id sp, int rep, int year);

	void createOccupancy(species_id sp, int nbOutputRows);

	bool outOccupancyHeaders(Species* pSpecies);
	void closeOccupancyOfs(species_id sp);
	void outOccupancy(Species* pSpecies);

	// Functions to handle input and output
	int readLandscape(
		int filenum,		// fileNum == 0 for (first) habitat file and optional patch file
							// fileNum > 0  for subsequent habitat files under the %cover option
		string habitatFileName,
		const map<species_id, string>& patchFileNames
	);
	int readCosts(const map<species_id, string>& costFiles);
	void resetVisits();
	void outVisits(species_id sp, int rep, int landNb);	// save SMS path visits map to raster text file
#if RS_RCPP
	void outPathsFinishReplicate();
	void outPathsStartReplicate(int);
#endif

private:
	bool isArtificial;
	bool usesPatches;
	bool spDist;			// initial species distribution loaded
	bool isFractal;	
	bool isContinuous;
	bool isDynamic;			// landscape changes during simulation
	bool habsAreIndexed;	// habitat codes have been converted to index numbers
	short rasterType;		// 0 = habitat codes 1 = % cover 2 = quality 9 = artificial landscape
	int landNum;			// landscape number
	int resol;				// cell size (m)
	int nHab;				// no. of habitats
	int nHabMax;			// max. no. of habitats (used for batch input only)
	int dimX, dimY;			// dimensions
	int minX, minY;			// minimum available X and Y co-ordinates (bottom left corner)
	int maxX, maxY;			// maximum available X and Y co-ordinates
	float minPct, maxPct;	// min and max percentage of habitat in a cell
	float propSuit;			// proportion of suitable cells
	float hurst;			// Hurst exponent
	int maxCells;			// max. cells per patch (artificial landscapes)
	double minEast;			// real world min co-ordinates
	double minNorth;		// read from habitat raster

	// list of cells in the landscape
	// cells MUST be loaded in the sequence ascending x within descending y
	Cell*** cells;

	// list of patches in the landscape - can be in any sequence
	map<species_id, vector<Patch*>> patchesList;

	// list of habitat codes
	std::vector<int> habCodes;

	// list of initial individual species distributions
	map<species_id, InitDist> distns;

	// list of cells to be initialised for ALL species
	std::vector<DistCell*> initcells;

	// patch connectivity matrices (one per species)
	// indexed by [start patch seq num][end patch seq num]
	map<species_id, int**> connectMatrices;
	map<species_id, ofstream> outConnMatrices;
	map<species_id, ofstream> outOccupOfs;

	// global environmental stochasticity (epsilon)
	float* epsGlobal;	// pointer to time-series	
	
	// list of dynamic landscape changes
	std::vector<landChange> landChanges;
	map<species_id, vector<patchChange>> patchChanges;
	map<species_id, vector<costChange>> costsChanges;
	// patch and costs change matrices (temporary - used when reading dynamic landscape)
	// indexed by [descending y][x]
	map<species_id, vector<vector<cellChange>>> patchChgMatrices;
	map<species_id, vector<vector<cellChange>>> costsChgMatrices;

};

// NOTE: the following function is not a behaviour of Landscape, as it is run by the
// batch routine to check raster files before any Landscape has been initiated
rasterdata CheckRasterFile(string);

void ReadSpDynLandFile(
	ifstream& ifsSpLand,
	map<species_id, string>& pathsToPatchMaps,
	map<species_id, string>& pathsToCostMaps,
	const int& nbSpecies
);

extern paramStoch* paramsStoch;
extern paramSim* paramsSim;
extern RSrandom* pRandom;

#ifdef UNIT_TESTS
landParams createDefaultLandParams(const int& dim);
void testLandscape();
#endif

#if RS_RCPP
extern rasterdata landRaster, patchraster, spdistraster, costsraster;
extern void EOFerrorR(string);
extern void StreamErrorR(string);
#endif

//---------------------------------------------------------------------------
#endif
