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

 RangeShifter v2.0 Patch

 Implements the class: Patch

 A patch is a collection of one or more Cells in the the gridded Landscape,
 which together provide the area in which a single demographic unit of a Species,
 i.e. a Population, can reproduce. One or more Populations (of different Species)
 form a Sub-community associated with the Patch.

 There is no requirement that all the Cells be adjacent, although in practice
 that would usually be the case.

 Each Patch must have a unique positive integer id number supplied by the user,
 and the matrix, i.e. any part of the landscape which is not a breeding patch,
 is represented by Patch 0. However, as patch numbers need not be sequential,
 an internal sequential number is also applied.

 For a 'cell-based model', the user supplies no patch numbers, and a separate
 Patch is generated internally for each Cell, i.e. the 'cell-based model' is a
 special case of the 'patch-based model' in which each Patch has a single Cell.
 Moreover, there is also the 'matrix' Patch 0, which has no cells, but which
 holds the disperser population whilst its Individuals are in transit.

 In a true patch-based model, each Patch hold a list of its constituent Cells,
 EXCEPT for the matrix Patch 0. This is because that list would be extremely
 long for a very large landscape in which suitable patches are small and/or rare,
 and removing Cells from it if the landscape is dynamic would be inefficient.

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 25 June 2021 by Steve Palmer

 ------------------------------------------------------------------------------*/

#ifndef PatchH
#define PatchH

#include <vector>
using namespace std;

#include "Parameters.h"
#include "Cell.h"
#include "Species.h"
#include "Patch.h"

class Population;
struct popStats;

 //---------------------------------------------------------------------------

struct patchLimits {
	int xMin, xMax, yMin, yMax;
};

class Patch {
public:
	Patch(int seqnum, int num);
	~Patch();

	int getSeqNum();
	int getPatchNum();
	int getNCells();

	patchLimits getLimits(); // Returns the minimum and maximum co-ordinates of the patch
	bool withinLimits( // Does the patch fall (partially) within a specified rectangle?
		patchLimits rect// structure holding the SW and NE co-ordinates of the rectangle
	);
	void resetLimits(); // Reset minimum and maximum co-ordinates of the patch
	
	void addCell(Cell* pCell, int x, int y);
	// Return co-ordinates of a specified cell
	locn getCellLocn(int ix);
	// Return pointer to a specified cell
	Cell* getCell(int ix);
	// Return co-ordinates of patch centroid
	locn getCentroid(); 
	void removeCell(Cell* pCell);
	Cell* getRandomCell();

	void setPop(Population* p);
	Population* getPop();
	void resetPop();

	// Record the presence of a potential settler within the Patch
	void incrPossSettler(Species* pSpecies, int sex);
	// Get number of a potential settlers within the Patch
	int getPossSettlers(Species* pSpecies, int sex);
	void resetPossSettlers();
	
	// Calculate total Patch carrying capacity (no. of inds)
	void setCarryingCapacity(Species* pSpecies, patchLimits landlimits,
		float epsGlobal, short nHab, short rasterType, short landIx, bool gradK);
	float getK();

	int getInitNbInds(const bool& isPatchModel, const int& landResol) const;
	
	float getEnvVal(const bool& isPatchModel, const float& epsGlobal);

	bool speciesIsPresent();

	void createOccupancy(int nbOutputRows);
	void updateOccupancy(int whichRow);
	int getOccupancy(int whichRow);

#ifndef NDEBUG
	// Testing only
	void overrideK(const float& k) { localK = k; }
#endif

private:
	int patchSeqNum;// sequential patch number - patch 0 is reserved for the inter-patch matrix
	int patchNum; 	// patch number as supplied by the user (not forced to be sequential)
	int nCells;			// no. of cells in the patch
	int xMin, xMax, yMin, yMax; 	// min and max cell co-ordinates
	int x, y;				// centroid co-ordinates (approx.)
	Population* pPop; // pointer to population associated with the patch
	float localK;		// patch carrying capacity (individuals)
	bool changed;
	short nTemp[gMaxNbSexes];	// no. of potential settlers in each sex
	vector<int> occupancy;	// pointer to occupancy array

	std::vector <Cell*> cells;
};

//---------------------------------------------------------------------------

extern paramStoch* paramsStoch;
extern RSrandom* pRandom;

#endif
