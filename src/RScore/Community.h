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
 It holds a SubCommunity for each Patch in the Landscape (including the matrix),
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

#include "SubCommunity.h"
#include "Landscape.h"
#include "Patch.h"
#include "Cell.h"
#include "Species.h"
#include "NeutralStatsManager.h"

//---------------------------------------------------------------------------
struct commStats {
	int ninds, nnonjuvs, suitable, occupied;
	int minX, maxX, minY, maxY;
};

class Community {

public:
	Community(Landscape*);
	~Community(void);
	SubCommunity* addPop(Patch*, int);
	// functions to manage populations occurring in the community
	void initialise(
		Species*,	// pointer to Species
		int				// year (relevent only for seedType == 2)
	);
	void addManuallySelected();
	void resetPopns();
	void initialInd(
		Landscape* pLandscape, 
		Species* pSpecies, 
		Patch* pPatch, 
		Cell* pCell, 
		int ix
	);
	void localExtinction(int option);
	void patchChanges();
	void reproduction(
		int				// year
	);
	void emigration(void);
	void dispersal(
		short landIx,
		short nextseason
	);

	// Remove emigrants from patch 0 (matrix) and transfer to SubCommunity in which
	// their destination co-ordinates fall (executed for the matrix patch only)
	void completeDispersal(Landscape* pLandscape, bool connect);

	void survival(
		short,	// part:		0 = determine survival & development,
		//		 			1 = apply survival changes to the population
		short,	// option0:	0 = stage 0 (juveniles) only         )
		//					1 = all stages                       ) used by part 0 only
		//					2 = stage 1 and above (all non-juvs) )
		short 	// option1:	0 - development only (when survival is annual)
		//		  	 		1 - development and survival
	);
	void ageIncrement();
	int totalInds();
	commStats getStats();

	void createOccupancy(int nbOutputRows, int nbReps);
	void updateOccupancy(int whichRow, int replicate);
	// Open occupancy file, write header record and set up occupancy array
	bool outOccupancyHeaders(
		int		// option: -999 to close the file
	);
	void outOccupancy();

	bool outRangeHeaders( // Open range file and write header record
		Species*,	// pointer to Species
		int				// Landscape number (-999 to close the file)
	);
	void outRange( // Write record to range file
		Species*, // pointer to Species
		int,			// replicate
		int,			// year
		int				// generation
	);

	bool outPopHeaders( // Open population file and write header record
		Species*, // pointer to Species
		int       // option: -999 to close the file
	);

	void outPop( // Write records to population file
		int,	// replicate
		int,	// year
		int		// generation
	);

	void outInds( // Write records to individuals file
		int,	// replicate
		int,	// year
		int,	// generation
		int		// Landscape number (>= 0 to open the file, -999 to close the file
					//									 -1 to write data records)
	);
	
	void outOccSuit(
		bool	// TRUE if occupancy graph is to be viewed on screen
	);
	bool outTraitsHeaders(
		Landscape* pLandscape,
		Species* pSpecies,
		int	landnb
	);
	bool outTraitsRowsHeaders( // Open trait rows file and write header record
		Species*, // pointer to Species
		int       // Landscape number (-999 to close the file)
	);
	void outTraits( // Write records to traits file
		Species*,		// pointer to Species
		int,				// replicate
		int,				// year
		int					// generation
	);
	void writeTraitsRows( // Write records to trait rows file
		Species*,	// pointer to Species
		int,			// replicate
		int,			// year
		int,			// generation
		int,			// row number (Y cell co-ordinate)
		traitsums	// structure holding sums of trait genes for dispersal (see Population.h)
	);
#if RS_RCPP && !R_CMD
	Rcpp::IntegerMatrix addYearToPopList(int, int);
#endif

	//sample individuals for genetics (or could be used for anything)
	void sampleIndividuals(Species* pSpecies);

	bool openOutGenesFile(const bool& isDiploid, const int landNr, const int rep);
	void outputGeneValues(const int& year, const int& gen, Species* pSpecies);

	//control neutral stat output
	void outNeutralGenetics(Species* pSpecies, int rep, int yr, int gen, bool outWeirCockerham, bool outWeirHill);

	//file openers
	bool openNeutralOutputFile(Species* pSpecies, const int landNr);
	bool openPerLocusFstFile(Species* pSpecies, Landscape* pLandscape, const int landNr, const int rep);
	bool openPairwiseFstFile(Species* pSpecies, Landscape* pLandscape, const int landNr, const int rep);

	//file writers
	void writeNeutralOutputFile(int rep, int yr, int gen, bool outWeirCockerham, bool outWeirHill);
	void writePerLocusFstatFile(Species* pSpecies, const int yr, const int gen, const  int nAlleles, const int nLoci, set<int> const& patchList);
	void writePairwiseFstFile(Species* pSpecies, const int yr, const int gen, const  int nAlleles, const int nLoci, set<int> const& patchList);

private:
	Landscape* pLandscape;
	int indIx;				// index used to apply initial individuals
	vector<vector <int>> occSuit;	// occupancy of suitable cells / patches

	Population* matrixPop;
	std::vector <Population*> popns;

	unique_ptr<NeutralStatsManager> pNeutralStatistics;
};

extern paramSim* paramsSim;
extern paramInit* paramsInit;


//---------------------------------------------------------------------------
#endif
