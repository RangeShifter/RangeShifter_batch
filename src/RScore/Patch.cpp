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


 //---------------------------------------------------------------------------

#include "Patch.h"
#include "Population.h"

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

Patch::Patch(int seqnum, int num, species_id whichSpecies)
{
	pPop = nullptr;
	sp = whichSpecies;
	patchSeqNum = seqnum; 
	patchNum = num; 
	nCells = 0;
	xMin = yMin = 999999999; 
	xMax = yMax = 0; 
	x = y = 0;
	localK = 0.0;
	for (int sex = 0; sex < gMaxNbSexes; sex++) {
		nTemp[sex] = 0;
	}
	changed = false;
	gradVal = 1.0; // i.e. no effect
}

Patch::~Patch() {
	cells.clear();
}

int Patch::getSeqNum() { return patchSeqNum; }

int Patch::getPatchNum() { return patchNum; }

int Patch::getNCells() { return nCells; }

patchLimits Patch::getLimits() {
	patchLimits p;
	p.xMin = xMin; 
	p.xMax = xMax; 
	p.yMin = yMin;
	p.yMax = yMax;
	return p;
}

// Does the patch fall (partially) within a specified rectangle?
bool Patch::withinLimits(patchLimits rect) {
	locn loc;
	if (xMin <= rect.xMax && xMax >= rect.xMin && yMin <= rect.yMax && yMax >= rect.yMin) {
		// patch is within the rectangle UNLESS it is irregular in shape and lies at a corner
		// of the rectangle
		if ((xMin >= rect.xMin && xMax <= rect.xMax)
			|| (yMin >= rect.yMin && yMax <= rect.yMax)) {
			// patch lies within or along an edge of the initialistaion rectangle
			return true;
		}
		else {
			// check for any cell of the patch lying within the rectangle
			int ncells = cells.size();
			for (int i = 0; i < ncells; i++) {
				loc = getCellLocn(i);
				if (loc.x >= rect.xMin && loc.x <= rect.xMax
					&& loc.y >= rect.yMin && loc.y <= rect.yMax) {
					// cell lies within the rectangle
					return true;
				}
			}
		}
	}
	return false;
}

// Reset minimum and maximum co-ordinates of the patch if it has been changed
void Patch::resetLimits() {
	if (changed) {
		// remove any deleted cells
		std::vector <Cell*> newcells; // for all retained and added cells
		int ncells = static_cast<int>(cells.size());
		for (int i = 0; i < ncells; i++) {
			if (cells[i] != nullptr) {
				newcells.push_back(cells[i]);
			}
		}
		cells.clear();
		cells = newcells;
		// reset patch limits
		locn loc;
		xMin = yMin = 999999999; 
		xMax = yMax = 0;
		ncells = static_cast<int>(cells.size());
		for (int i = 0; i < ncells; i++) {
			loc = getCellLocn(i);
			if (loc.x < xMin) xMin = loc.x;
			if (loc.x > xMax) xMax = loc.x;
			if (loc.y < yMin) yMin = loc.y;
			if (loc.y > yMax) yMax = loc.y;
		}
		changed = false;
	}
}

// Add a cell to the patch
void Patch::addCell(Cell* pCell, int x, int y) {
	cells.push_back(pCell);
	nCells++;
	if (x < xMin) xMin = x;
	if (x > xMax) xMax = x;
	if (y < yMin) yMin = y;
	if (y > yMax) yMax = y;
}

// Calculate the total carrying capacity (no. of individuals) and
// centroid co-ordinates of the patch
void Patch::setCarryingCapacity(Species* pSpecies, float epsGlobal,
	short nHab, short rasterType, short dynLandIndex) {
	
	envStochParams env = paramsStoch->getStoch();
	locn loc;
	short hx;
	float k, q;

	localK = 0.0; // no. of suitable cells (unadjusted K > 0) in the patch
	int nsuitable = 0;
	double mean;

	if (!pSpecies->isWithinLimits(xMin, yMin) 
		&& !pSpecies->isWithinLimits(xMax, yMax)) {
		// patch lies wholly outwith current landscape limit
		return; // K is 0
	}

	int ncells = static_cast<int>(cells.size());
	int xsum = 0, ysum = 0;

	// Calculate environmental contribution 
	double envVal = gradVal; // 0 if no gradient
	if (env.usesStoch && env.inK) {
		envVal += env.stochIsLocal ?
			cells[0]->getEps() : // cell-based if stoch is local so only one cell/patch
			epsGlobal;
	}

	for (int i = 0; i < ncells; i++) {

		switch (rasterType) {
		case 0: // habitat codes
			hx = cells[i]->getHabIndex(dynLandIndex);
			k = pSpecies->getHabK(hx);
			if (k > 0.0) {
				nsuitable++;
				localK += envVal * k;
			}
			break;
		case 1: // cover %
			k = 0.0;
			for (int j = 0; j < nHab; j++) { // loop through cover layers
				q = cells[i]->getHabitat(j);
				k += q * pSpecies->getHabK(j) / 100.0f;
			}
			if (k > 0.0) {
				nsuitable++;
				localK += envVal * k;
			}
			break;
		case 2: // habitat quality
			q = cells[i]->getHabitat(dynLandIndex);
			if (q > 0.0) {
				nsuitable++;
				localK += envVal * pSpecies->getHabK(0) * q / 100.0f;
			}
			break;
		}
		loc = cells[i]->getLocn();
		xsum += loc.x; 
		ysum += loc.y;
	}
// calculate centroid co-ordinates
	if (ncells > 0) {
		mean = static_cast<double>(xsum) / ncells;
		x = static_cast<int>(mean + 0.5);
		mean = static_cast<double>(ysum) / ncells;
		y = static_cast<int>(mean + 0.5);
	}
	if (env.usesStoch && env.inK) { // environmental stochasticity in K
		// apply min and max limits to K over the whole patch
		// NB limits have been stored as N/cell rather than N/ha
		float limit;
		limit = pSpecies->getMinMax(0) * nsuitable;
		if (localK < limit) localK = limit;
		limit = pSpecies->getMinMax(1) * nsuitable;
		if (localK > limit) localK = limit;
	}
}

float Patch::getK() { return localK; }

bool Patch::isSuitable() { return localK > 0.0f; }

int Patch::getInitNbInds(initParams init, const bool& isPatchModel, const int& landResol) const {

	int nInds = 0;
	if (localK > 0.0) { // patch is currently suitable for this species
		switch (init.initDens) {
		case 0: // at carrying capacity
			nInds = static_cast<int>(localK);
			break;
		case 1: // at half carrying capacity
			nInds = static_cast<int>(localK / 2.0);
			break;
		case 2: // specified no. per cell or density
			if (isPatchModel) {
				nInds = static_cast<int>(init.indsHa * (float)(nCells * landResol * landResol) / 10000.0);
			}
			else {
				nInds = init.indsCell * nCells;
			}
			break;
		}
	}
	return nInds;
}

void Patch::calcGradVal(Species* pSpecies) {
	Cell* pCell = getCell(0); // only one cell/patch since cell-based landscape
	envGradParams grad = pSpecies->getEnvGradient();
	float dev = pCell->getEnvDev();
	float distFromOpt = fabs(grad.optY - y);
	gradVal = min(1.0, 1.0 - distFromOpt * grad.gradIncr + dev * grad.factor);
	if (gradVal < 0.000001) gradVal = 0.0;
}

double Patch::getGradVal() const {
	return gradVal;
}


// Return co-ordinates of a specified cell
locn Patch::getCellLocn(int ix) {
	locn loc; 
	loc.x = -666; 
	loc.y = -666;
	if (ix >= 0 && ix < cells.size()) {
		loc = cells[ix]->getLocn();
	}
	return loc;
}
// Return pointer to a specified cell
Cell* Patch::getCell(int ix) {
	if (ix >= 0 && ix < cells.size()) return cells[ix];
	else return nullptr;
}
// Return co-ordinates of patch centroid
locn Patch::getCentroid() {
	locn loc; 
	loc.x = x; 
	loc.y = y;
	return loc;
}

// Select a Cell within the Patch at random, and return pointer to it
// For a cell-based model, this will be the only Cell
Cell* Patch::getRandomCell() {
	Cell* pCell = nullptr;
	int ix;
	int ncells = static_cast<int>(cells.size());
	if (ncells > 0) {
		if (ncells == 1) ix = 0;
		else ix = pRandom->IRandom(0, ncells - 1);
		pCell = cells[ix];
	}
	return pCell;
}

// Remove a cell from the patch
void Patch::removeCell(Cell* pCell) {
	int ncells = (int)cells.size();
	for (int i = 0; i < ncells; i++) {
		if (pCell == cells[i]) {
			cells[i] = nullptr; 
			i = ncells;
			nCells--;
			changed = true;
		}
	}
}

void Patch::setPop(Population* p) {
	pPop = p;
}

Population* Patch::getPop() {
	return pPop;
}

void Patch::resetPop() {
	pPop = nullptr;
}

void Patch::resetPossSettlers() {
	for (int sex = 0; sex < gMaxNbSexes; sex++) {
		nTemp[sex] = 0;
	}
}

// Record the presence of a potential settler within the Patch
void Patch::incrPossSettler(int sex) {
	nTemp[sex]++;
}

// Get number of a potential settlers within the Patch
int Patch::getPossSettlers(int sex) {
	return nTemp[sex];
}

bool Patch::speciesIsPresent() {
	return getPop() != nullptr;
}

void Patch::createOccupancy(int nbOutputRows) {
	occupancy = vector<int>(nbOutputRows, 0);
}

void Patch::updateOccupancy(int whichRow) {
	popStats ps = pPop->getStats();
	if (ps.nInds > 0 && ps.breeding)
		occupancy[whichRow]++;
}

int Patch::getOccupancy(int whichRow) {
	return occupancy[whichRow];
}

// Increment the tally of overlapping cells with the patch of another species
void Patch::incrementPatchOverlap(Patch* pOverlappingPatch) {

	// Add entries if not already present
	const species_id otherSp = pOverlappingPatch->getSpeciesID();
	const int otherPatchId = pOverlappingPatch->getPatchNum();
	if (!overlappingPatches.contains(otherSp))
		overlappingPatches.emplace(otherSp, map<int, double>());
	auto& thisSpEntry = overlappingPatches.at(sp);
	if (!thisSpEntry.contains(otherPatchId))
		thisSpEntry.emplace(otherPatchId, 0.0f);

	// Increment overlap (in nb of cells) between these patches;
	thisSpEntry.at(otherPatchId)++;
}

// Calculate final amount of overlap with other species patches
// as proportion of this patch's cells that overlap
void Patch::calcPatchOverlap() {

	for (auto& patchesFromThisSp : overlappingPatches) {
		for (auto& thisOtherPatch : patchesFromThisSp.second)
			// convert the previously calculated tally into a proportion of overlap
			thisOtherPatch.second /= nCells;
	}

}

void Patch::resetPatchOverlap() {
	overlappingPatches.clear();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------



