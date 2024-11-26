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

 RangeShifter v2.0 Community

 Implements the Community class

 There is ONLY ONE instance of a Community in an individual replicate simulation.
 It holds a Population for each Patch in the Landscape (including the matrix),
 and is thus the highest-level entity accessed for most processing concerned with
 simulated populations.

 Optionally, the Community maintains a record of the occupancy of suitable cells
 or patches during the course of simulation of multiple replicates.

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 25 June 2021 by Anne-Kathleen Malchow

 ------------------------------------------------------------------------------*/

#ifndef CommunityH
#define CommunityH

#include <vector>
#include <algorithm>
#include <memory>
#include <ranges>
using namespace std;

#include "Landscape.h"
#include "Patch.h"
#include "Cell.h"
#include "Species.h"
#include "NeutralStatsManager.h"
#include "Parameters.h"
#include "Population.h"

//---------------------------------------------------------------------------
struct commStats {
	int ninds, nnonjuvs, suitable, occupied;
	int minX, maxX, minY, maxY;
};

class Community {

public:
	Community(Landscape* pLand, map<int, Species*> allSpecies);
	~Community();
	// functions to manage populations occurring in the community
	void initialise(Species* pSpecies, int year);
	void resetPopns();
	Species* findSpecies(int speciesID);
	void initialInd(
		Landscape* pLandscape, 
		Species* pSpecies, 
		Patch* pPatch, 
		Cell* pCell, 
		int ix
	);
	void localExtinction(int option);
	void patchChanges();
	void reproduction(int year);
	void emigration();
	void dispersal(short landIx, short nextseason);

	// Remove emigrants from patch 0 (matrix) and transfer to the Population in which
	// their destination co-ordinates fall (executed for the matrix patch only)
	void completeDispersal(Landscape* pLandscape, bool connect);

	void drawSurvivalDevlpt(
		bool resolveJuvs,
		bool resolveAdults,
		bool resolveDev,
		bool resolveSurv
	);
	void applySurvivalDevlpt();

	void ageIncrement();
	int totalInds();
	commStats getStats();

	void createOccupancy(int nbOutputRows, int nbReps);
	void updateOccupancy(int whichRow, int replicate);

	// Open occupancy file, write header record and set up occupancy array
	bool outOccupancyHeaders();
	void outOccupancy();
	bool closeOccupancyOfs();
	void outOccSuit(bool view);

	bool outRangeHeaders( // Open range file and write header record
		Species* pSpecies,
		int	landnr // (-999 to close the file)
	);
	// Write record to range file
	void outRange(Species* pSpecies, int rep, int yr, int gen);
	bool closeRangeOfs();

	// Open population file and write header record
	bool outPopHeaders(Species* pSpecies);
	bool closePopOfs();

	// Write records to population file
	void outPop(int rep, int year, int gen);

	void outIndsHeaders(int rep, int landNr, bool patchModel, Species* pSpecies);
	void closeOutIndsOfs();

	// Write records to individuals file
	void outInds(int rep, int year,	int gen);
	
	// Write records to traits file
	void outTraits(Species* pSpecies, int rep, int year, int gen);
	bool outTraitsHeaders(Landscape* pLandscape, Species* pSpecies, int landnb);
	bool closeOutTraitOfs();

	// Open trait rows file and write header record
	bool outTraitsRowsHeaders(Species* pSpecies, int landnr);
	// Write records to trait rows file
	void writeTraitsRows(
		Species* pSpecies, 
		int rep, 
		int year,
		int gen, 
		int row, 
		traitsums ts
	);
	bool closeTraitRows();

#if RS_RCPP && !R_CMD
	Rcpp::IntegerMatrix addYearToPopList(int rep, int yr);
#endif

	// sample individuals for genetics (or could be used for anything)
	void sampleIndividuals(Species* pSpecies);

	bool openOutGenesFile(const bool& isDiploid, const int landNr, const int rep);
	void outputGeneValues(const int& year, const int& gen, Species* pSpecies);
	bool closeOutGenesOfs();

	// control neutral stat output
	void outNeutralGenetics(Species* pSpecies, int rep, int yr, int gen, bool outWeirCockerham, bool outWeirHill);
	bool openNeutralOutputFile(Species* pSpecies, const int landNr);
	void writeNeutralOutputFile(int rep, int yr, int gen, bool outWeirCockerham, bool outWeirHill);
	bool closeNeutralOutputOfs();

	bool openPerLocusFstFile(Species* pSpecies, Landscape* pLandscape, const int landNr, const int rep);
	bool closePerLocusFstFile();
	void writePerLocusFstatFile(Species* pSpecies, const int yr, const int gen, const  int nAlleles, const int nLoci, set<int> const& patchList);

	bool openPairwiseFstFile(Species* pSpecies, Landscape* pLandscape, const int landNr, const int rep);
	bool closePairwiseFstFile();
	void writePairwiseFstFile(Species* pSpecies, const int yr, const int gen, const  int nAlleles, const int nLoci, set<int> const& patchList);

private:
	map<int, Species*> speciesMap;
	Landscape* pLandscape;
	int indIx;				// index used to apply initial individuals
	vector<vector <int>> occSuit;	// occupancy of suitable cells / patches

	Population* matrixPop;
	std::vector <Population*> popns;

	unique_ptr<NeutralStatsManager> pNeutralStatistics;

	ofstream outPopOfs;
	ofstream outIndsOfs;
	ofstream outOccupOfs;
	ofstream outSuitOfs;
	ofstream outTraitsOfs;
	ofstream outRangeOfs;
	ofstream outTraitsRows;
	ofstream ofsGenes;
	ofstream outPairwiseFstOfs;
	ofstream outWCFstatOfs;
	ofstream outPerLocusFstat;
};

extern paramSim* paramsSim;
extern paramInit* paramsInit;


//---------------------------------------------------------------------------
#endif
