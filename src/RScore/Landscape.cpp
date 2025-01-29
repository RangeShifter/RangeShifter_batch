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

#include "Landscape.h"
//---------------------------------------------------------------------------

ifstream landscape;

#if RS_RCPP
ofstream outMovePaths;
#endif // RS_RCPP

//---------------------------------------------------------------------------

bool isInLandscape(const int& x, const int& y, const landData& land) {
	return (x >= land.minX
		&& x <= land.maxX
		&& y >= land.minY
		&& y <= land.maxY);
}

// Initial species distribution functions

InitDist::InitDist(Species* pSp)
{
	pSpecies = pSp;
	resol = 1;
	maxX = 0;
	maxY = 0;
	minEast = 0.0;
	minNorth = 0.0;
}

InitDist::~InitDist() {
	int ncells = cells.size();
	for (int i = 0; i < ncells; i++)
		if (cells[i] != nullptr) delete cells[i];
	cells.clear();
}

void InitDist::setDistribution(int nInit) {
	int rr = 0;
	int ncells = (int)cells.size();
	if (nInit == 0) { // set all cells to be initialised
		for (int i = 0; i < ncells; i++) {
			cells[i]->setCell(true);
		}
	}
	else { // set specified number of cells at random to be initialised
		if (nInit > ncells / 2) { // use backwards selection method
			for (int i = 0; i < ncells; i++) cells[i]->setCell(true);
			for (int i = 0; i < (ncells - nInit); i++) {
				do {
					rr = pRandom->IRandom(0, ncells - 1);
				} while (!cells[rr]->selected());
				cells[rr]->setCell(false);
			}
		}
		else { // use forwards selection method
			for (int i = 0; i < ncells; i++) cells[i]->setCell(false);
			for (int i = 0; i < nInit; i++) {
				do {
					rr = pRandom->IRandom(0, ncells - 1);
				} while (cells[rr]->selected());
				cells[rr]->setCell(true);
			}
		}
	}
}

// Set a specified cell (by position in cells vector)
void InitDist::setDistCell(int ix, bool init) {
	cells[ix]->setCell(init);
}

// Set a specified cell (by co-ordinates)
void InitDist::setDistCell(locn loc, bool init) {
	locn cellloc;
	int ncells = (int)cells.size();
	for (int i = 0; i < ncells; i++) {
		cellloc = cells[i]->getLocn();
		if (cellloc.x == loc.x && cellloc.y == loc.y) {
			cells[i]->setCell(init);
			i = ncells;
		}
	}
}

// Specified location is within the initial distribution?
bool InitDist::isInInitialDist(locn loc) {
	int ncells = static_cast<int>(cells.size());
	for (int i = 0; i < ncells; i++) {
		if (cells[i]->toInitialise(loc)) { // cell is to be initialised
			return true;
		}
	}
	return false;
}

int InitDist::cellCount() {
	return cells.size();
}

// Return the co-ordinates of a specified initial distribution cell
locn InitDist::getCell(int ix) {
	locn loc;
	if (ix >= 0 && ix < cells.size()) {
		loc = cells[ix]->getLocn();
	}
	else {
		loc.x = loc.y = -666; // indicates invalid index specified
	}
	return loc;
}

// Return the co-ordinates of a specified initial distribution cell if it has been
// selected - otherwise return negative co-ordinates
locn InitDist::getSelectedCell(int ix) {
	locn loc; 
	loc.x = loc.y = -666;
	if (ix < cells.size()) {
		if (cells[ix]->selected()) {
			loc = cells[ix]->getLocn();
		}
	}
	return loc;
}

locn InitDist::getDimensions() {
	locn d; 
	d.x = maxX; 
	d.y = maxY; 
	return d;
}

void InitDist::resetDistribution() {
	int ncells = static_cast<int>(cells.size());
	for (int i = 0; i < ncells; i++) {
		cells[i]->setCell(false);
	}
}

//---------------------------------------------------------------------------

// Read species initial distribution file

#if RS_RCPP
int InitDist::readDistribution(string distfile) {
	wstring header;
	int patchCode, nodata;
	int ncols, nrows;
	wifstream dfile; // species distribution file input stream

	// open distribution file
	dfile.open(distfile, std::ios::binary);
	if (spdistraster.utf) {
		// apply BOM-sensitive UTF-16 facet
		dfile.imbue(std::locale(dfile.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
	}
	if (!dfile.is_open()) return 21;

	// read landscape data from header records of distribution file
	// NB headers of all files have already been compared
	dfile >> header >> ncols
		>> header >> nrows
		>> header >> minEast
		>> header >> minNorth
		>> header >> resol
		>> header >> nodata;

	if (!dfile.good()) {
		// corrupt file stream
		StreamErrorR(distfile);
		dfile.close();
		dfile.clear();
		return 144;
	}

	maxX = ncols - 1;
	maxY = nrows - 1;

	// set up bad integer value to ensure that valid values are read
	int badvalue = -9;
	if (nodata == -9) badvalue = -99;

	for (int y = nrows - 1; y >= 0; y--) {
		for (int x = 0; x < ncols; x++) {
			patchCode = badvalue;
			if (dfile >> patchCode) {
				if (patchCode == 1) { // species present
					cells.push_back(new DistCell(x, y));
				}
				else if (patchCode != nodata && patchCode != 0) { // error in file
					dfile.close();
					dfile.clear();
					return 22;
				}
			}
			else {
				// corrupt file stream
#if !R_CMD
				Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
				StreamErrorR(distfile);
				dfile.close();
				dfile.clear();
				return 144;
			}
		}
	}
	dfile >> patchCode;
	if (!dfile.eof()) EOFerrorR(distfile);

	dfile.close();
	dfile.clear();
	return 0;
}

#else

int InitDist::readDistribution(string distfile) {
	string header;
	int p, nodata;
	int ncols, nrows;
	ifstream dfile; // species distribution file input stream

	// open distribution file
	dfile.open(distfile.c_str());
	if (!dfile.is_open()) return 21;

	// read landscape data from header records of distribution file
	// NB headers of all files have already been compared
	dfile >> header >> ncols
		>> header >> nrows
		>> header >> minEast
		>> header >> minNorth
		>> header >> resol
		>> header >> nodata;

	maxX = ncols - 1;
	maxY = nrows - 1;

	// set up bad integer value to ensure that valid values are read
	int badvalue = -9;
	if (nodata == -9) badvalue = -99;

	for (int y = nrows - 1; y >= 0; y--) {
		for (int x = 0; x < ncols; x++) {
			p = badvalue;
			dfile >> p;

			if (p == 1) { // species present
				cells.push_back(new DistCell(x, y));
			}
			else if (p == nodata || p == 0) {
				// nothing, species absent
			}
			else { // any otehr value is an error
				dfile.close();
				dfile.clear();
				return 22;
			}
		}
	}
	dfile.close();
	dfile.clear();
	return 0;
}

#endif // not RS_RCPP

//---------------------------------------------------------------------------

// Landscape functions

Landscape::Landscape(const speciesMap_t& allSpecies) {
	usesPatches = false; 
	spDist = false; 
	generated = false;
	fractal = false; 
	continuous = false;
	dynamic = false; 
	habsAreIndexed = false;
	resol = spResol = 1;
	landNum = 0;
	rasterType = 0;
	nHab = nHabMax = 0;
	dimX = dimY = 100;
	minX = minY = 0;
	maxX = maxY = 99;
	minPct = maxPct = propSuit = hurst = 0.0;
	maxCells = 100;
	minEast = minNorth = 0.0;
	cells = nullptr;
	epsGlobal = nullptr;

	// Initialise maps for species-dependent members
	for (auto& sp : views::keys(allSpecies)) {
		patchesList.emplace(sp, vector<Patch*>());
		patchChgMatrices.emplace(sp, nullptr);
		costsChgMatrices.emplace(sp, nullptr);
		connectMatrices.emplace(sp, nullptr);
		patchChanges.emplace(sp, vector<patchChange>());
		costsChanges.emplace(sp, vector<costChange>());
	}
}

Landscape::~Landscape() {

	if (cells != nullptr) {
		for (int y = dimY - 1; y >= 0; y--) {
			for (int x = 0; x < dimX; x++) {
				if (cells[y][x] != nullptr) delete cells[y][x];
			}
			if (cells[y] != nullptr) {
				delete[] cells[y];
			}
		}
		delete[] cells;
		cells = nullptr;
	}

	for (auto& [sp, patches] : patchesList) {
		int npatches = static_cast<int>(patches.size());
		for (int i = 0; i < npatches; i++)
			if (patches[i] != nullptr) delete patches[i];
		patches.clear();
	}
	
	int ndistns = static_cast<int>(distns.size());
	for (int i = 0; i < ndistns; i++)
		if (distns[i] != nullptr) delete distns[i];
	distns.clear();

	int ninitcells = static_cast<int>(initcells.size());
	for (int i = 0; i < ninitcells; i++)
		if (initcells[i] != nullptr) delete initcells[i];
	initcells.clear();

	habCodes.clear();
	landChanges.clear();
	for (const species_id& sp : views::keys(patchChanges))
		patchChanges.at(sp).clear();
	for (const species_id& sp : views::keys(connectMatrices))
		deleteConnectMatrix(sp);
	for (const species_id& sp : views::keys(patchChgMatrices))
		deletePatchChgMatrix(sp);
	if (epsGlobal != nullptr) delete[] epsGlobal;
}

// Remove all patches and cells
// Used for replicating artificial landscape without deleting the landscape itself
void Landscape::resetLand() {
	resetLandLimits();

	for (auto& [sp, patches] : patchesList) {
		int npatches = static_cast<int>(patches.size());
		for (int i = 0; i < npatches; i++)
			if (patches[i] != nullptr) delete patches[i];
		patches.clear();
	}
	
	if (cells != nullptr) {
		for (int y = dimY - 1; y >= 0; y--) {
			for (int x = 0; x < dimX; x++) {
				if (cells[y][x] != nullptr) delete cells[y][x];
			}
			if (cells[y] != nullptr) {
				delete[] cells[y];
			}
		}
		delete[] cells;
		cells = nullptr;
	}
}

void Landscape::setLandParams(landParams ppp, bool batchMode)
{
	generated = ppp.generated; 
	usesPatches = ppp.usesPatches; 
	spDist = ppp.useSpDist;
	dynamic = ppp.dynamic;
	landNum = ppp.landNum;
	if (ppp.resol > 0) resol = ppp.resol;
	if (ppp.spResol > 0 
		&& ppp.spResol % ppp.resol == 0) spResol = ppp.spResol;
	if ((ppp.rasterType >= 0 
		&& ppp.rasterType <= 2) || ppp.rasterType == 9)
		rasterType = ppp.rasterType;
	if (ppp.nHab >= 1) nHab = ppp.nHab;
	if (ppp.nHabMax >= 1) nHabMax = ppp.nHabMax;
	if (ppp.dimX > 0) dimX = ppp.dimX;
	if (ppp.dimY > 0) dimY = ppp.dimY;
	if (ppp.minX >= 0 
		&& ppp.maxX >= 0 
		&& ppp.minX <= ppp.maxX 
		&& ppp.maxX < dimX) {
		minX = ppp.minX; 
		maxX = ppp.maxX;
	}
	else {
		minX = 0; 
		maxX = dimX - 1;
	}
	if (ppp.minY >= 0 
		&& ppp.maxY >= 0 
		&& ppp.minY <= ppp.maxY
		&& ppp.maxY < dimY) {
		minY = ppp.minY; 
		maxY = ppp.maxY;
	}
	else {
		minY = 0; 
		maxY = dimY - 1;
	}
	if (batchMode && rasterType == 0) {
		// in batch mode, set up sequential habitat codes if not already present
		if (habCodes.size() == 0) {
			for (int i = 0; i < nHabMax; i++) {
				habCodes.push_back(i + 1);
			}
		}
	}
}

landParams Landscape::getLandParams() {
	landParams ppp;
	ppp.generated = generated; 
	ppp.usesPatches = usesPatches; 
	ppp.useSpDist = spDist;
	ppp.dynamic = dynamic;
	ppp.landNum = landNum;
	ppp.resol = resol;
	ppp.spResol = spResol;
	ppp.rasterType = rasterType;
	ppp.nHab = nHab; 
	ppp.nHabMax = nHabMax;
	ppp.dimX = dimX; 
	ppp.dimY = dimY;
	ppp.minX = minX; 
	ppp.minY = minY;
	ppp.maxX = maxX; 
	ppp.maxY = maxY;
	return ppp;
}

landData Landscape::getLandData() {
	landData dd;
	dd.resol = resol;
	dd.dimX = dimX; 
	dd.dimY = dimY;
	dd.minX = minX; 
	dd.minY = minY;
	dd.maxX = maxX;
	dd.maxY = maxY;
	return dd;
}

void Landscape::setGenLandParams(genLandParams ppp)
{
	fractal = ppp.fractal;
	continuous = ppp.continuous;
	if (ppp.minPct > 0.0 && ppp.minPct < 100.0) minPct = ppp.minPct;
	if (ppp.maxPct > 0.0 && ppp.maxPct <= 100.0) maxPct = ppp.maxPct;
	if (ppp.propSuit >= 0.0 && ppp.propSuit <= 1.0) propSuit = ppp.propSuit;
	if (ppp.hurst > 0.0 && ppp.hurst < 1.0) hurst = ppp.hurst;
	if (ppp.maxCells > 0) maxCells = ppp.maxCells;
}

genLandParams Landscape::getGenLandParams()
{
	genLandParams ppp;
	ppp.fractal = fractal;
	ppp.continuous = continuous;
	ppp.minPct = minPct; 
	ppp.maxPct = maxPct; 
	ppp.propSuit = propSuit; 
	ppp.hurst = hurst;
	ppp.maxCells = maxCells;
	return ppp;
}

void Landscape::setLandLimits(int x0, int y0, int x1, int y1) {
	if (x0 >= 0 && x1 >= 0 && x0 <= x1 && x1 < dimX
		&& y0 >= 0 && y1 >= 0 && y0 <= y1 && y1 < dimY) {
		minX = x0; maxX = x1; minY = y0; maxY = y1;
	}
}

void Landscape::resetLandLimits() {
	minX = minY = 0; 
	maxX = dimX - 1; 
	maxY = dimY - 1;
}

//---------------------------------------------------------------------------

landOrigin Landscape::getOrigin() {
	landOrigin originVal;
	originVal.minEast = minEast; 
	originVal.minNorth = minNorth;
	return originVal;
}

//---------------------------------------------------------------------------

// Functions to handle habitat codes

bool Landscape::habitatsIndexed() { return habsAreIndexed; }

void Landscape::listHabCodes() {
	int nhab = static_cast<int>(habCodes.size());
#if RS_RCPP && !R_CMD
	Rcpp::Rcout << endl;
	for (int i = 0; i < nhab; i++) {
		Rcpp::Rcout << "Habitat code[ " << i << "] = " << habCodes[i] << endl;
	}
	Rcpp::Rcout << endl;
#else
	cout << endl;
	for (int i = 0; i < nhab; i++) {
		cout << "Habitat code[ " << i << "] = " << habCodes[i] << endl;
	}
	cout << endl;
#endif
}

void Landscape::addHabCode(int hab) {
	int nhab = static_cast<int>(habCodes.size());
	for (int i = 0; i < nhab; i++) {
		if (hab == habCodes[i]) return; // already exists
	}
	habCodes.push_back(hab); 
	nHab++;
}

// Get the index number of the specified habitat in the habitats vector
int Landscape::findHabCode(int hab) {
	int nhab = static_cast<int>(habCodes.size());
	for (int i = 0; i < nhab; i++) {
		if (hab == habCodes[i]) return i;
	}
	return -999;
}

// Get the specified habitat code
int Landscape::getHabCode(int ixhab) {
	if (ixhab < static_cast<int>(habCodes.size())) return habCodes[ixhab];
	else return -999;
}

void Landscape::clearHabitats() {
	habCodes.clear();
}

//---------------------------------------------------------------------------
void Landscape::setCellArray() {
	if (cells != nullptr) resetLand();
	cells = new Cell* *[dimY];
	for (int y = dimY - 1; y >= 0; y--) {
		cells[y] = new Cell* [dimX];
		for (int x = 0; x < dimX; x++) {
			cells[y][x] = nullptr;
		}
	}
}

//---------------------------------------------------------------------------
/* Create an artificial landscape (random or fractal), which can be
either binary (habitat index 0 is the matrix, 1 is suitable habitat)
or continuous (0 is the matrix, >0 is suitable habitat) */
void Landscape::generatePatches(const speciesMap_t& allSpecies)
{
	int x, y, ncells;
	double p;
	Patch* pPatch;
	Cell* pCell;

	setCellArray();
	// as landscape generator returns cells in a random sequence, first set up all cells
	// in the landscape in the correct sequence, then update them and create patches for
	// habitat cells
	for (int yy = dimY - 1; yy >= 0; yy--) {
		for (int xx = 0; xx < dimX; xx++) {
			addNewCellToLand(xx, yy, 0);
		}
	}

	rasterType = continuous ? rasterType = 2 : 0;

	int patchnum = 1;

	for (auto& sp : views::keys(allSpecies)) {

		vector<Patch*> patches;
		// Each species has a matrix patch with index 0
		Patch* matrixPatch = new Patch(0, 0, sp);
		patches.push_back(matrixPatch);

		if (fractal) {
			p = 1.0 - propSuit;
			// createFractalLandscape() requires Max_prop > 1 (but does not check it!)
			// as in turn it calls runif(1.0,Max_prop)
			double maxpct = maxPct < 1.0 ? 100.0 : maxPct;

			vector<fractalPatch> FracLandscape = createFractalLandscape(dimY, dimX, hurst, p, maxpct, minPct);

			vector<fractalPatch>::iterator iter = FracLandscape.begin();
			while (iter != FracLandscape.end()) {
				x = iter->y_coord;
				y = iter->x_coord;
				pCell = findCell(x, y);
				if (continuous) {
					if (iter->value > 0.0) { // habitat
						pPatch = new Patch(patchnum, patchnum, sp);
						patchnum++;
						patches.push_back(pPatch);
						addCellToPatch(sp, pCell, pPatch, iter->value);
					}
					else { // matrix
						addCellToPatch(sp, pCell, matrixPatch, iter->value);
					}
				}
				else { // discrete
					if (iter->avail == 0) { // matrix
						addCellToPatch(sp, pCell, matrixPatch);
					}
					else { // habitat
						pPatch = new Patch(patchnum, patchnum, sp);
						patchnum++;
						patches.push_back(pPatch);
						addCellToPatch(sp, pCell, pPatch);
						pCell->changeHabIndex(0, 1);
					}
				}
				iter++;
			}
		}
		else { // random landscape
			int hab = 0;
			ncells = static_cast<int>(dimX * dimY * propSuit + 0.00001); // no. of cells to initialise
			int i = 0;
			do {
				do {
					x = pRandom->IRandom(0, dimX - 1);
					y = pRandom->IRandom(0, dimY - 1);
					pCell = findCell(x, y);
					hab = pCell->getHabIndex(0);
				} while (hab > 0);

				pPatch = new Patch(patchnum, patchnum, sp);
				patchnum++;
				patches.push_back(pPatch);
				addCellToPatch(sp, findCell(x, y), pPatch);
				pCell->changeHabIndex(0, 1);

				if (continuous) {
					pCell->addHabitat((float)(minPct + pRandom->Random() * (maxPct - minPct)));
				}
				i++;
			} while (i < ncells);

			// remaining cells need to be added to the matrix patch
			p = 0.0;
			x = 0;
			for (int yy = dimY - 1; yy >= 0; yy--) {
				for (int xx = 0; xx < dimX; xx++) {

					pCell = findCell(xx, yy);
					if (continuous && pCell->getHabitat(0) <= 0.0) {
						addCellToPatch(sp, pCell, matrixPatch, (float)p);
					}
					else if (pCell->getHabIndex(0) == 0) {
						addCellToPatch(sp, pCell, matrixPatch, x);
					}

				}
			}

		} // fractal or not

		patchesList.at(sp) = patches;

	} // loop through species
}

//---------------------------------------------------------------------------

// Landscape patch-management functions

//---------------------------------------------------------------------------
/* Create a patch for each suitable cell of a cell-based landscape (all other
habitat cells are added to the matrix patch) */
void Landscape::allocatePatches(const speciesMap_t& allSpecies)
{
	float habK;
	Patch* pPatch;
	Cell* pCell;

	// Delete all existing patches (from previous landscape)
	for (auto& [sp, patches] : patchesList) {
		for (int i = 0; i < patches.size(); i++) {
			if (patches[i] != nullptr) delete patches[i];
		}
		patches.clear();
	}
	
	int patchnum = 1; // patch number is unique across all species

	for (auto& [sp, pSpecies] : allSpecies) {

		vector<Patch*> patches;

		// Create one matrix patch for each species
		Patch* matrixPatch = new Patch(0, 0, sp); // all matrix patches share index 0
		patches.push_back(matrixPatch);

		switch (rasterType) {

		case 0: // habitat codes

			for (int y = dimY - 1; y >= 0; y--) {
				for (int x = 0; x < dimX; x++) {

					if (cells[y][x] != nullptr) { // not no-data cell
						pCell = cells[y][x];
						habK = 0.0;
						int nhab = pCell->nHabitats();
						for (int i = 0; i < nhab; i++) {
							habK += pSpecies->getHabK(pCell->getHabIndex(i));
						}
						if (habK > 0.0) { // cell is suitable - create a patch for it
							pPatch = new Patch(patchnum, patchnum, sp);
							patchnum++;
							patches.push_back(pPatch);
							addCellToPatch(sp, pCell, pPatch);
						}
						else { // cell is not suitable - add to the matrix patch
							addCellToPatch(sp, pCell, matrixPatch);
						}
					}
				}
			}
			break;
		case 1: // habitat cover

			for (int y = dimY - 1; y >= 0; y--) {
				for (int x = 0; x < dimX; x++) {

					if (cells[y][x] != nullptr) { // not no-data cell
						pCell = cells[y][x];
						habK = 0.0;
						int nhab = pCell->nHabitats();
						for (int i = 0; i < nhab; i++) {
							habK += pSpecies->getHabK(i) * pCell->getHabitat(i) / 100.0f;
						}
						if (habK > 0.0) { // cell is suitable - create a patch for it
							pPatch = new Patch(patchnum, patchnum, sp);
							patchnum++;
							patches.push_back(pPatch);
							addCellToPatch(sp, pCell, pPatch);
						}
						else { // cell is not suitable - add to the matrix patch
							addCellToPatch(sp, pCell, matrixPatch);
						}
					}
				}
			}
			break;
		case 2: // habitat quality
			for (int y = dimY - 1; y >= 0; y--) {
				for (int x = 0; x < dimX; x++) {
					if (cells[y][x] != nullptr) { // not no-data cell
						pCell = cells[y][x];
						habK = 0.0;
						int nhab = pCell->nHabitats();
						for (int i = 0; i < nhab; i++) {
							habK += pSpecies->getHabK(0) * pCell->getHabitat(i) / 100.0f;
						}
						if (habK > 0.0) { // cell is suitable (at some time) - create a patch for it
							pPatch = new Patch(patchnum, patchnum, sp);
							patchnum++;
							patches.push_back(pPatch);
							addCellToPatch(sp, pCell, pPatch);
						}
						else { // cell is never suitable - add to the matrix patch
							addCellToPatch(sp, pCell, matrixPatch);
						}
					}
				}
			}
			break;

		} // end of switch (rasterType)

		patchesList.at(sp) = patches;

	} // end of loop through species
}

Patch* Landscape::addNewPatch(species_id id, int num, species_id sp)
{
	patchesList.at(id).push_back(new Patch(num, num, sp));
	return patchesList.at(id)[patchesList.at(id).size() - 1];
}

Patch* Landscape::addNewPatch(species_id id, int seqnum, int num, species_id sp)
{
	patchesList.at(id).push_back(new Patch(seqnum, num, sp));
	return patchesList.at(id)[patchesList.at(id).size() - 1];
}

void Landscape::resetPatchLimits() {
	for (auto& [sp, patches] : patchesList) {
		int npatches = static_cast<int>(patches.size());
		for (int i = 0; i < npatches; i++) {
			patches[i]->resetLimits();
		}
	}
}

void Landscape::addNewCellToLand(int x, int y, float q) {
	if (q < 0.0) // no-data cell - no Cell created
		cells[y][x] = nullptr;
	else {
		set<species_id> spLabels;
		for (const species_id& sp : views::keys(patchesList)) {
			spLabels.insert(sp);
		}
		cells[y][x] = new Cell(x, y, nullptr, q, spLabels);
	}
}

void Landscape::addNewCellToLand(int x, int y, int hab) {
	if (hab < 0) // no-data cell - no Cell created
		cells[y][x] = nullptr;
	else {
		set<species_id> spLabels;
		for (const species_id& sp : views::keys(patchesList)) {
			spLabels.insert(sp);
		}
		cells[y][x] = new Cell(x, y, nullptr, hab, spLabels);
	}
}

void Landscape::addCellToLand(Cell* c) {
	if (cells == 0) throw runtime_error("Landscape cells member is uninitialised.");
	if (c->getHabIndex(0) < 0.0)
		throw logic_error("Can't add no-data cell to landscape.");
	locn l = c->getLocn();
	cells[l.y][l.x] = c;
}

void Landscape::addNewCellToPatch(Patch* pPatch, int x, int y, float q) {
	if (q < 0.0) throw logic_error("Attempt to add a cell with negative habitat quality.");
	else { // create the new cell
		set<species_id> spLabels;
		for (const species_id& sp : views::keys(patchesList)) {
			spLabels.insert(sp);
		}
		cells[y][x] = new Cell(x, y, pPatch, q, spLabels);
		if (pPatch != nullptr) { // not the matrix patch
			// add the cell to the patch
			pPatch->addCell(cells[y][x], x, y);
		}
	}
}

void Landscape::addNewCellToPatch(Patch* pPatch, int x, int y, int hab) {
	if (hab < 0) throw logic_error("Attempt to add a cell with negative habitat code.");
	else { // create the new cell
		set<species_id> spLabels;
		for (const species_id& sp : views::keys(patchesList)) {
			spLabels.insert(sp);
		}
		cells[y][x] = new Cell(x, y, pPatch, hab, spLabels);
		if (pPatch != nullptr) { // not the matrix patch
			// add the cell to the patch
			pPatch->addCell(cells[y][x], x, y);
		}
	}
}

void Landscape::addCellToPatch(species_id whichSpecies, Cell* pCell, Patch* pPatch) {
	pCell->setPatch(whichSpecies, pPatch);
	locn loc = pCell->getLocn();
	// add the cell to the patch
	pPatch->addCell(pCell, loc.x, loc.y);
}

void Landscape::addCellToPatch(species_id whichSpecies, Cell* pCell, Patch* pPatch, float q) {
	pCell->setPatch(whichSpecies, pPatch);
	// update the habitat type of the cell
	pCell->addHabitat(q);
	locn loc = pCell->getLocn();
	// add the cell to the patch
	pPatch->addCell(pCell, loc.x, loc.y);
}

void Landscape::addCellToPatch(species_id whichSpecies, Cell* pCell, Patch* pPatch, int hab) {
	pCell->setPatch(whichSpecies, pPatch);
	// update the habitat type of the cell
	pCell->addHabIndex(hab);
	locn loc = pCell->getLocn();
	// add the cell to the patch
	pPatch->addCell(pCell, loc.x, loc.y);
}

patchData Landscape::getPatchData(species_id id, int patchID) {
	patchData ppp;
	Patch* pPatch = patchesList.at(id)[patchID];
	ppp.pPatch = pPatch;
	ppp.patchNum = pPatch->getPatchNum();
	ppp.nCells = pPatch->getNCells();
	locn randloc;
	randloc.x = -666; 
	randloc.y = -666;
	Cell* pCell = pPatch->getRandomCell();
	if (pCell != nullptr) {
		randloc = pCell->getLocn();
	}
	ppp.x = randloc.x; 
	ppp.y = randloc.y;
	return ppp;
}

Patch* Landscape::findPatch(species_id whichSpecies, int patchID) {
	for (auto p : patchesList.at(whichSpecies)) {
		if (patchID == p->getPatchNum()) return p;
	};
	return nullptr;
}

bool Landscape::existsPatch(species_id whichSpecies, int patchID) {
	return findPatch(whichSpecies, patchID) != nullptr;
}

void Landscape::samplePatches(speciesMap_t& allSpecies, const string& samplingOption) {

	for (auto& [sp, pSpecies] : allSpecies) {

		vector<int> sampledPatches;
		vector<int> eligiblePatches;
		int nbToSample = pSpecies->getNbPatchesToSample();

		// Get list of viable patches where the species is present
		for (auto p : patchesList.at(sp)) {
			if (p->isMatrix()) continue; // skip
			if (samplingOption == "random") { // then all patches are eligible
				eligiblePatches.push_back(p->getPatchNum());
			}
			else if (p->speciesIsPresent()) {
				// only patches with at least 1 ind can be sampled
				eligiblePatches.push_back(p->getPatchNum());
			}
		}

		if (samplingOption == "all") {
			sampledPatches = eligiblePatches;
		}
		else if (samplingOption == "random_occupied" || samplingOption == "random") {
			if (nbToSample > eligiblePatches.size())
				nbToSample = eligiblePatches.size();
			auto rng = pRandom->getRNG();
			sample(eligiblePatches.begin(), eligiblePatches.end(), std::back_inserter(sampledPatches),
				nbToSample, rng);
		}
		else {
			throw logic_error("Sampling option should be random, rnadom_occupied or all when sampling patches.");
		}

		set<int> patchIds;
		copy(sampledPatches.begin(), sampledPatches.end(), inserter(patchIds, patchIds.end()));
		pSpecies->setSamplePatchList(patchIds);

	} // end loop through species
}

void Landscape::resetPatchPopns() {
	for (auto& [sp, patches] : patchesList) {
		int npatches = static_cast<int>(patches.size());
		for (int i = 0; i < npatches; i++) {
			patches[i]->resetPop();
		}
	}
}

void Landscape::updateCarryingCapacity(const speciesMap_t& allSpecies, int yr, short landIx) {
	
	envGradParams grad = paramsGrad->getGradient();
	bool gradK = false;
	if (grad.gradient && grad.gradType == 1) gradK = true; // gradient in carrying capacity
	patchLimits landlimits;
	landlimits.xMin = minX; 
	landlimits.xMax = maxX;
	landlimits.yMin = minY; 
	landlimits.yMax = maxY;

	for (auto& [sp, patches] : patchesList) {
		int npatches = static_cast<int>(patches.size());
		for (int i = 0; i < npatches; i++) {
			if (!patches[i]->isMatrix()) {
				patches[i]->setCarryingCapacity(allSpecies.at(sp), landlimits, getGlobalStoch(yr),
					nHab, rasterType, landIx, gradK);
			}
		}
	}
}

Cell* Landscape::findCell(int x, int y) {
	if (x >= 0 && x < dimX && y >= 0 && y < dimY) return cells[y][x];
	else return 0;
}

int Landscape::patchCount(species_id id) const {
	return static_cast<int>(patchesList.at(id).size());
}

int Landscape::allPatchCount() const {
	int count = 0;
	for (const species_id& sp : views::keys(patchesList))
		count += patchCount(sp);
	return count;
}

// Check that total cover of any cell does not exceed 100%
// and identify matrix cells
int Landscape::checkTotalCover() {
	if (rasterType != 1) return 0; // not appropriate test
	int nCells = 0;
	for (int y = dimY - 1; y >= 0; y--) {
		for (int x = 0; x < dimX; x++) {
			if (cells[y][x] != nullptr)
			{ // not a no-data cell
				float sumCover = 0.0;
				for (int i = 0; i < nHab; i++) {
					sumCover += cells[y][x]->getHabitat(i);
				}
				if (sumCover > 100.00001) nCells++; // decimal part to allow for floating point error
				if (sumCover <= 0.0) // cell is a matrix cell
					cells[y][x]->addHabIndex(0);
				else
					cells[y][x]->addHabIndex(1);
			}
		}
	}
	return nCells;
}

// Convert habitat codes stored on loading habitat codes landscape to
// sequential sorted index numbers
void Landscape::updateHabitatIndices() {
	// sort codes
	sort(habCodes.begin(), habCodes.end());
	nHab = static_cast<int>(habCodes.size());
	// convert codes in landscape
	int habIx;
	int nbLandChanges = static_cast<int>(landChanges.size());
	for (int y = dimY - 1; y >= 0; y--) {
		for (int x = 0; x < dimX; x++) {

			if (cells[y][x] != nullptr) { // not a no-data cell
				for (int c = 0; c <= nbLandChanges; c++) {
					habIx = cells[y][x]->getHabIndex(c);
					if (habIx >= 0) {
						habIx = findHabCode(habIx);
						cells[y][x]->changeHabIndex(c, habIx);
					}
				}
			}
		}
	}
	habsAreIndexed = true;
}

void Landscape::setEnvGradient(Species* pSpecies, bool initial)
{
	float dist_from_opt, dev;
	float habK;
	//int hab;
	double envval;
	// gradient parameters
	envGradParams grad = paramsGrad->getGradient();
	for (int y = dimY - 1; y >= 0; y--) {
		for (int x = 0; x < dimX; x++) {
			// NB: gradient lies in range 0-1 for all types, and is applied when necessary...
			// ... implies gradient increment will be dimensionless in range 0-1 (but << 1)
			if (cells[y][x] != 0) { // not no-data cell
				habK = 0.0;
				int nhab = cells[y][x]->nHabitats();
				for (int i = 0; i < nhab; i++) {
					switch (rasterType) {
					case 0:
						habK += pSpecies->getHabK(cells[y][x]->getHabIndex(i));
						break;
					case 1:
						habK += pSpecies->getHabK(i) * cells[y][x]->getHabitat(i) / 100.0f;
						break;
					case 2:
						habK += pSpecies->getHabK(0) * cells[y][x]->getHabitat(i) / 100.0f;
						break;
					}
				}
				if (habK > 0.0) { // suitable cell
					if (initial) { // set local environmental deviation
						cells[y][x]->setEnvDev((float)pRandom->Random() * (2.0f) - 1.0f);
					}
					dist_from_opt = (float)(fabs((double)grad.opt_y - (double)y));
					dev = cells[y][x]->getEnvDev();
					envval = 1.0 - dist_from_opt * grad.grad_inc + dev * grad.factor;
					if (envval < 0.000001) envval = 0.0;
					if (envval > 1.0) envval = 1.0;
				}
				else envval = 0.0;
				cells[y][x]->setEnvVal((float)envval);
			}
		}
	}

}

void Landscape::setGlobalStoch(int nyears) {
	envStochParams env = paramsStoch->getStoch();
	if (epsGlobal != 0) delete[] epsGlobal;
	epsGlobal = new float[nyears];
	epsGlobal[0] = (float)(pRandom->Normal(0.0, env.std) * sqrt(1.0 - (env.ac * env.ac)));
	for (int i = 1; i < nyears; i++) {
		epsGlobal[i] = (float)(env.ac * epsGlobal[i - 1] + pRandom->Normal(0.0, env.std) * sqrt(1.0 - (env.ac * env.ac)));
	}
}

float Landscape::getGlobalStoch(int yr) {
	if (epsGlobal != 0 && yr >= 0) {
		return epsGlobal[yr];
	}
	else return 0.0;
}

void Landscape::updateLocalStoch(void) {
	envStochParams env = paramsStoch->getStoch();
	float randpart;
	for (int y = dimY - 1; y >= 0; y--) {
		for (int x = 0; x < dimX; x++) {
			if (cells[y][x] != 0) { // not a no-data cell
				randpart = (float)(pRandom->Normal(0.0, env.std) * sqrt(1.0 - (env.ac * env.ac)));
				cells[y][x]->updateEps((float)env.ac, randpart);
			}
		}
	}

}

void Landscape::resetCosts(void) {
	for (int y = dimY - 1; y >= 0; y--) {
		for (int x = 0; x < dimX; x++) {
			if (cells[y][x] != 0) { // not a no-data cell
				cells[y][x]->resetCost();
			}
		}
	}
}

void Landscape::resetEffCosts(void) {
	for (int y = dimY - 1; y >= 0; y--) {
		for (int x = 0; x < dimX; x++) {
			if (cells[y][x] != 0) { // not a no-data cell
				cells[y][x]->resetEffCosts();
			}
		}
	}
}

//---------------------------------------------------------------------------

// Dynamic landscape functions

void Landscape::setDynamicLand(bool dyn) { dynamic = dyn; }

void Landscape::addLandChange(landChange c) {
	landChanges.push_back(c);
}

int Landscape::numLandChanges(void) { return (int)landChanges.size(); }

landChange Landscape::getLandChange(short ix) {
	landChange c; c.chgnum = c.chgyear = 0;
	c.habfile = c.pchfile = c.costfile = "none";
	int nchanges = (int)landChanges.size();
	if (ix < nchanges) c = landChanges[ix];
	return c;
}

void Landscape::deleteLandChanges(void) {
	while (landChanges.size() > 0) landChanges.pop_back();
	landChanges.clear();
}

#if RS_RCPP && !R_CMD
int Landscape::readLandChange(int filenum, bool usesCosts, wifstream& ifsHabMap, wifstream& ifsPatchMap, wifstream& ifsDynCostFile, int noDataHabCode, int noDataPatch, int costNoData) {
#else
int Landscape::readLandChange(int filenum, bool usesCosts) {
#endif

#if RS_RCPP
	wstring header;
#else
	string header;
	int ncols, nrows, habNoData, costNoData, patchNoData;
	costNoData = 0;
	patchNoData = 0;
#endif
	int habCode = 0, patchCode = 0, costCode = 0;
	float habFloat, patchFloat, costFloat;
	simParams sim = paramsSim->getSim();

	if (filenum < 0) return 19;
	const species_id sp = gSingleSpeciesID; // only one species for now
	int patchSeq = usesPatches ? patchCount(sp) : 0;

#if !RS_RCPP || R_CMD
	ifstream ifsDynHabFile; // habitat file input stream
	ifstream ifsDynPatchFile; // patch file input stream
	ifstream ifsDynCostFile; // costs file input stream

	// open habitat file and optionally also patch and costs files
	ifsDynHabFile.open(landChanges[filenum].habfile.c_str());
	if (!ifsDynHabFile.is_open()) return 30;
	if (usesPatches) {
		ifsDynPatchFile.open(landChanges[filenum].pchfile.c_str());
		if (!ifsDynPatchFile.is_open()) {
			ifsDynHabFile.close(); 
			ifsDynHabFile.clear();
			return 31;
		}
	}
	if (usesCosts) {
		ifsDynCostFile.open(landChanges[filenum].costfile.c_str());
		if (!ifsDynCostFile.is_open()) {
			ifsDynHabFile.close(); ifsDynHabFile.clear();
			if (ifsDynPatchFile.is_open()) {
				ifsDynPatchFile.close(); ifsDynPatchFile.clear();
			}
			return 32;
		}
	}

	// read header records of habitat (and patch) file(s)
	// NB headers of all files have already been compared
	ifsDynHabFile
		>> header >> ncols 
		>> header >> nrows 
		>> header >> habFloat 
		>> header >> habFloat
		>> header >> habFloat 
		>> header >> habNoData;
	if (usesPatches) {
		for (int i = 0; i < 5; i++) ifsDynPatchFile 
			>> header >> patchFloat;
		ifsDynPatchFile 
			>> header >> patchNoData;
	}
	if (usesCosts) {
		for (int i = 0; i < 5; i++) ifsDynCostFile 
			>> header >> costFloat;
		ifsDynCostFile 
			>> header >> costNoData;
	}
#endif

	// set up bad float values to ensure that valid values are read
	float badHabFloat = -9.0; 
	if (habNoData == -9) badHabFloat = -99.0;
	float badPatchFloat = -9.0;
	if (patchNoData == -9) badPatchFloat = -99.0;
	float badCostFloat = -9.0; 
	if (costNoData == -9) badCostFloat = -99.0;

	switch (rasterType) {

	case 0: // raster with habitat codes - 100% habitat each cell
		for (int y = dimY - 1; y >= 0; y--) {
			for (int x = 0; x < dimX; x++) {
				habFloat = badHabFloat;
#if RS_RCPP
				if (ifsHabMap >> habFloat) {
					habCode = static_cast<int>(habFloat);
				}
				else {
					// corrupt file stream
#if !R_CMD
					Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
					StreamErrorR("habitatchgfile");
					ifsHabMap.close();
					ifsHabMap.clear();
					ifsPatchMap.close();
					ifsPatchMap.clear();
					return 171;
				}

				if (usesPatches) {
					patchFloat = badPatchFloat;
					if (ifsPatchMap >> patchFloat) {
						patchCode = static_cast<int>(patchFloat);
					}
					else {
						// corrupt file stream
#if !R_CMD
						Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
						StreamErrorR("patchchgfile");
						ifsHabMap.close();
						ifsHabMap.clear();
						ifsPatchMap.close();
						ifsPatchMap.clear();
						return 172;
					}

				}
				if (usesCosts) {
					costFloat = badCostFloat;
					if (ifsDynCostFile >> costFloat) {
						costCode = (int)costFloat;
					}
					else {
						// corrupt file stream
#if !R_CMD
						Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
						StreamErrorR("costchgfile");
						ifsHabMap.close();
						ifsHabMap.clear();
						ifsPatchMap.close();
						ifsPatchMap.clear();
						return 173;
					}
				}
#else // not RCPP
				ifsDynHabFile >> habFloat;
				habCode = static_cast<int>(habFloat);
				if (usesPatches) {
					patchFloat = badPatchFloat;
					ifsDynPatchFile >> patchFloat;
					patchCode = static_cast<int>(patchFloat);

					}
				if (usesCosts) {
					costFloat = badCostFloat;
					ifsDynCostFile >> costFloat;
					costCode = static_cast<int>(costFloat);

				}
#endif
				if (cells[y][x] != nullptr) { // not a no data cell (in initial landscape)
					if (habCode == habNoData) { // invalid no data cell in change map
						ifsDynHabFile.close();
						ifsDynHabFile.clear();
						return 36;
					}
					else if (habCode < 0 
						|| (sim.batchMode && (habCode < 1 || habCode > nHabMax))) {
						// invalid habitat code
						ifsDynHabFile.close();
						ifsDynHabFile.clear();
						if (usesPatches) {
							ifsDynPatchFile.close();
							ifsDynPatchFile.clear();
						}
						return 33;
					}
					else { // valid habitat code
						addHabCode(habCode);
						cells[y][x]->addHabIndex(habCode);
					}

					if (usesPatches) {
						if (patchCode < 0 || patchCode == patchNoData) { // invalid patch code
#if RS_RCPP && !R_CMD
							if (patchCode == noDataPatch) Rcpp::Rcout << "Found patch NA in valid habitat cell." << std::endl;
							else Rcpp::Rcout << "Found negative patch ID in valid habitat cell." << std::endl;
#endif
							ifsDynHabFile.close();
							ifsDynHabFile.clear();
							ifsDynPatchFile.close();
							ifsDynPatchFile.clear();
							return 34;
						}
						else {
							patchChgMatrices.at(sp)[y][x].nextVal = patchCode;
							if (patchCode > 0 && !existsPatch(sp, patchCode)) {
								patchesList.at(sp).push_back(new Patch(patchSeq, patchCode, sp));
								patchSeq++;
							}
						}
					}
					if (usesCosts) {
						if (costCode < 1) { // invalid cost
#if RS_RCPP
							Rcpp::Rcout << "Found invalid cost value of " << changeIndex << " in cell x " << x << " and y  " << y << std::endl;
#endif
							ifsDynHabFile.close();
							ifsDynHabFile.clear();
							if (ifsDynPatchFile.is_open()) {
								ifsDynPatchFile.close();
								ifsDynPatchFile.clear();
							}
							return 38;
						}
						else {
							costsChgMatrices.at(sp)[y][x].nextVal = costCode;
						}
					}
				} // if cell exists

			} // for x
		} // for y

#if RS_RCPP
		ifsHabMap >> habFloat;
		if (!ifsHabMap.eof()) EOFerrorR("habitatchgfile");
		if (usesPatches) {
			ifsPatchMap >> patchFloat;
			if (!ifsPatchMap.eof()) EOFerrorR("patchchgfile");
		}
		if (usesCosts) {
			ifsDynCostFile >> costFloat;
			if (!ifsDynCostFile.eof()) EOFerrorR("costchgfile");
		}
#endif
		break;

	case 2: // habitat quality
		for (int y = dimY - 1; y >= 0; y--) {
			for (int x = 0; x < dimX; x++) {
				habFloat = badHabFloat;

#if RS_RCPP
				if (ifsHabMap >> habFloat) {
					habCode = static_cast<int>(habFloat);
				}
				else {
					// corrupt file stream
#if !R_CMD
					Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
					StreamErrorR("habitatchgfile");
					ifsHabMap.close();
					ifsHabMap.clear();
					ifsPatchMap.close();
					ifsPatchMap.clear();
					return 172;
				}
				if (usesPatches) {
					patchFloat = badPatchFloat;
					if (ifsPatchMap >> patchFloat) {
						patchCode = static_cast<int>(patchFloat);
					}
					else {
						// corrupt file stream
#if !R_CMD
						Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
						StreamErrorR("patchchgfile");
						ifsHabMap.close();
						ifsHabMap.clear();
						ifsPatchMap.close();
						ifsPatchMap.clear();
						return 175;
					}
				}
				if (usesCosts) {
					costFloat = badCostFloat;
					if (ifsDynCostFile >> costFloat) {
						costCode = (int)costFloat;
					}
					else {
						// corrupt file stream
#if RS_RCPP && !R_CMD
						Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
						StreamErrorR("costchgfile");
						ifsHabMap.close();
						ifsHabMap.clear();
						ifsPatchMap.close();
						ifsPatchMap.clear();
						return 173;
					}
				}
#else
				ifsDynHabFile >> habFloat;
				habCode = static_cast<int>(habFloat);
				if (usesPatches) {
					patchFloat = badPatchFloat;
					ifsDynPatchFile >> patchFloat;
					patchCode = static_cast<int>(patchFloat);
				}
				if (usesCosts) {
					costFloat = badCostFloat;
					ifsDynCostFile >> costFloat;
					costCode = static_cast<int>(costFloat);
				}
#endif
				if (cells[y][x] != nullptr) { // not a no data cell (in initial landscape)
					if (habCode == habNoData) { // invalid no data cell in change map
						ifsDynHabFile.close();
						ifsDynHabFile.clear();
						if (usesPatches) {
							ifsDynPatchFile.close();
							ifsDynPatchFile.clear();
						}
						return 36;
					}
					else {
						if (habFloat < 0.0 || habFloat > 100.0) { // invalid quality score
							ifsDynHabFile.close();
							ifsDynHabFile.clear();
							if (usesPatches) {
								ifsDynPatchFile.close();
								ifsDynPatchFile.clear();
							}
							return 37;
						}
						else {
							cells[y][x]->addHabitat(habFloat);
						}
					}
					if (usesPatches) {
						if (patchCode < 0 || patchCode == patchNoData) { // invalid patch code
#if RS_RCPP && !R_CMD
							if (patchCode == noDataPatch) Rcpp::Rcout << "Found patch NA in valid habitat cell." << std::endl;
							else Rcpp::Rcout << "Found negative patch ID in valid habitat cell." << std::endl;
#endif
							ifsDynHabFile.close();
							ifsDynHabFile.clear();
							ifsDynPatchFile.close();
							ifsDynPatchFile.clear();
							return 34;
						}
						else {
							patchChgMatrices.at(sp)[y][x].nextVal = patchCode;
							if (patchCode > 0 && !existsPatch(sp, patchCode)) {
								// Create the patch if it doesn't exist already
								patchesList.at(sp).push_back(new Patch(patchSeq, patchCode, sp));
								patchSeq++;
							}
						}
					}
					if (usesCosts) {
						if (costCode < 1) { // invalid cost
#if RS_RCPP
							Rcpp::Rcout << "Found invalid cost value of " << changeIndex << "in cell x " << x << " and y  " << y << std::endl;
#endif
							ifsDynHabFile.close(); 
							ifsDynHabFile.clear();
							if (ifsDynPatchFile.is_open()) {
								ifsDynPatchFile.close(); 
								ifsDynPatchFile.clear();
							}
							return 38;
						}
						else {
							costsChgMatrices.at(sp)[y][x].nextVal = costCode;
						}
					}
				}
			} // end x
		} // end y
#if RS_RCPP
		ifsHabMap >> habFloat;
		if (!ifsHabMap.eof()) EOFerrorR("habitatchgfile");
		if (usesPatches) {
			ifsPatchMap >> patchFloat;
			if (!ifsPatchMap.eof()) EOFerrorR("patchchgfile");
		}
		if (usesCosts) {
			ifsDynCostFile >> costFloat;
			if (!ifsDynCostFile.eof()) EOFerrorR("costchgfile");
		}
#endif
		break;

	default:
		break;
	}

	if (ifsDynHabFile.is_open()) { ifsDynHabFile.close(); ifsDynHabFile.clear(); }
	if (ifsDynPatchFile.is_open()) { ifsDynPatchFile.close(); ifsDynPatchFile.clear(); }
	if (ifsDynCostFile.is_open()) { ifsDynCostFile.close(); ifsDynCostFile.clear(); }
	return 0;

}

// Create & initialise patch change matrix
void Landscape::createPatchChgMatrix() {
	Patch* pPatch;
	Cell* pCell;

	for (auto& [sp, patchChangeMatrix] : patchChgMatrices) {

		if (patchChangeMatrix != nullptr) deletePatchChgMatrix(sp);
		patchChangeMatrix = new cellChange* [dimY];

		for (int y = dimY - 1; y >= 0; y--) {

			patchChangeMatrix[y] = new cellChange[dimX];

			for (int x = 0; x < dimX; x++) {

				patchChangeMatrix[y][x] = cellChange();
				pCell = findCell(x, y);

				if (pCell == nullptr) { // no-data cell
					patchChangeMatrix[y][x].originVal = patchChangeMatrix[y][x].currentVal = 0;
				}
				else {
					// record initial patch number
					pPatch = pCell->getPatch(sp);
					if (pPatch == nullptr) { // matrix cell
						patchChangeMatrix[y][x].originVal = patchChangeMatrix[y][x].currentVal = 0;
					}
					else {
						patchChangeMatrix[y][x].originVal = patchChangeMatrix[y][x].currentVal = pPatch->getPatchNum();
					}
				}
				patchChangeMatrix[y][x].nextVal = 0;
			}
		}
	}
}

void Landscape::recordPatchChanges(int landIx) {

	patchChange chg;

	for (auto& [sp, patchChangeMatrix] : patchChgMatrices) {

		for (int y = dimY - 1; y >= 0; y--) {
			for (int x = 0; x < dimX; x++) {

				if (landIx == 0) { // reset to original landscape
					if (patchChangeMatrix[y][x].originVal != patchChangeMatrix[y][x].nextVal) {
						// record change of patch for current cell
						chg.chgnum = 666666;
						chg.x = x;
						chg.y = y;
						chg.oldpatch = patchChangeMatrix[y][x].nextVal;
						chg.newpatch = patchChangeMatrix[y][x].originVal;
						patchChanges.at(sp).push_back(chg);
					}
				}
				else { // any other change
					if (patchChangeMatrix[y][x].nextVal != patchChangeMatrix[y][x].currentVal) {
						// record change of patch for current cell
						chg.chgnum = landIx;
						chg.x = x;
						chg.y = y;
						chg.oldpatch = patchChangeMatrix[y][x].currentVal;
						chg.newpatch = patchChangeMatrix[y][x].nextVal;
						patchChanges.at(sp).push_back(chg);
					}
				}
				// reset cell for next landscape change
				patchChangeMatrix[y][x].currentVal = patchChangeMatrix[y][x].nextVal;
			}
		}
	}
}

int Landscape::numPatchChanges(species_id sp) { return static_cast<int>(patchChanges.at(sp).size()); }

patchChange Landscape::getPatchChange(species_id sp, int i) {
	return patchChanges.at(sp)[i];
}

void Landscape::deletePatchChgMatrices() {

	for (auto& sp : views::keys(patchChgMatrices)) {
		deletePatchChgMatrix(sp);
	}
}

void Landscape::deletePatchChgMatrix(species_id sp) {
	if (patchChgMatrices.at(sp) != nullptr) {
		for (int y = dimY - 1; y >= 0; y--) {
			delete[] patchChgMatrices.at(sp)[y];
		}
	}
	patchChgMatrices.at(sp) = nullptr;
}

// Create & initialise costs change matrix
void Landscape::createCostsChgMatrix()
{
	Cell* pCell;
	for (auto& [sp, costChangeMatrix] : costsChgMatrices) {
		if (costChangeMatrix != nullptr) deleteCostsChgMatrix(sp);
		costChangeMatrix = new cellChange * [dimY];
		for (int y = dimY - 1; y >= 0; y--) {
			costChangeMatrix[y] = new cellChange[dimX];
			for (int x = 0; x < dimX; x++) {
				costChangeMatrix[y][x] = cellChange();
				pCell = findCell(x, y);
				if (pCell == nullptr) { // no-data cell
					costChangeMatrix[y][x].originVal = costChangeMatrix[y][x].currentVal = 0;
				}
				else {
					// record initial cost
					costChangeMatrix[y][x].originVal = costChangeMatrix[y][x].currentVal = pCell->getCost();
				}
				costChangeMatrix[y][x].nextVal = 0;
			}
		}
	}
}

void Landscape::recordCostChanges(int landIx) {

	for (auto& [sp, costChangeMatrix] : costsChgMatrices) {

		if (costChangeMatrix == nullptr) return; // should not occur
		costChange chg;

		for (int y = dimY - 1; y >= 0; y--) {
			for (int x = 0; x < dimX; x++) {

				if (landIx == 0) { // reset to original landscape
					if (costChangeMatrix[y][x].originVal != costChangeMatrix[y][x].nextVal) {
						// record change of cost for current cell
						chg.chgnum = 666666;
						chg.x = x;
						chg.y = y;
						chg.oldcost = costChangeMatrix[y][x].nextVal;
						chg.newcost = costChangeMatrix[y][x].originVal;
						costsChanges.at(sp).push_back(chg);
					}
				}
				else { // any other change
					if (costChangeMatrix[y][x].nextVal != costChangeMatrix[y][x].currentVal) {
						// record change of cost for current cell
						chg.chgnum = landIx; chg.x = x; chg.y = y;
						chg.oldcost = costChangeMatrix[y][x].currentVal;
						chg.newcost = costChangeMatrix[y][x].nextVal;
						costsChanges.at(sp).push_back(chg);
					}
				}
				// reset cell for next landscape change
				costChangeMatrix[y][x].currentVal = costChangeMatrix[y][x].nextVal;

			} // for x
		} // for y 
	} // for species
}

int Landscape::getNbCostChanges(species_id sp) { return static_cast<int>(costsChanges.at(sp).size()); }

costChange Landscape::getCostChange(species_id sp, int i) {
	return costsChanges.at(sp)[i];
}

void Landscape::deleteCostsChgMatrices() {

	for (auto& sp : views::keys(costsChgMatrices)) {
		deleteCostsChgMatrix(sp);
	}
}

void Landscape::deleteCostsChgMatrix(species_id sp) {
	if (costsChgMatrices.at(sp) != nullptr) {
		for (int y = dimY - 1; y >= 0; y--) {
			delete[] costsChgMatrices.at(sp)[y];
		}
	}
	costsChgMatrices.at(sp) = nullptr;
}

//---------------------------------------------------------------------------

// Species distribution functions

int Landscape::newDistribution(Species* pSpecies, string distname) {
	int ndistns = distns.size();
	distns.push_back(new InitDist(pSpecies));
	int readcode = distns[ndistns]->readDistribution(distname);
	if (readcode != 0) { // error encountered
		// delete the distribution created above
		delete distns[ndistns];
		distns.pop_back();
	}
	return readcode;
}

void Landscape::setDistribution(Species* pSpecies, int nInit) {
	// WILL NEED TO SELECT DISTRIBUTION FOR CORRECT SPECIES ...
	// ... CURRENTLY IT IS THE ONLY ONE
	distns[0]->setDistribution(nInit);
}

// Specified cell match one of the distribution cells to be initialised?
bool Landscape::isInInitialDist(Species* pSpecies, locn loc) {
	// convert landscape co-ordinates to distribution co-ordinates
	locn initloc;
	initloc.x = loc.x * resol / spResol;
	initloc.y = loc.y * resol / spResol;
	// WILL HAVE TO GET CORRECT SPECIES WHEN THERE ARE MULTIPLE SPECIES ...
	bool initialise = distns[0]->isInInitialDist(initloc);
	return initialise;
}

void Landscape::deleteDistribution(Species* pSpecies) {
	if (distns[0] != nullptr) delete distns[0]; 
	distns.clear();
}

// Return no. of initial distributions
int Landscape::distnCount(void) {
	return (int)distns.size();
}

int Landscape::distCellCount(int dist) {
	return distns[dist]->cellCount();
}

// Set a cell in a specified initial distribution (by position in cells vector)
void Landscape::setDistnCell(int dist, int ix, bool init) {
	distns[dist]->setDistCell(ix, init);
}

// Set a cell in a specified initial distribution (by given co-ordinates)
void Landscape::setDistnCell(int dist, locn loc, bool init) {
	distns[dist]->setDistCell(loc, init);
}

// Get the co-ordinates of a specified cell in a specified initial distribution
locn Landscape::getDistnCell(int dist, int ix) {
	return distns[dist]->getCell(ix);
}

// Get the co-ordinates of a specified cell in a specified initial distribution
// Returns negative co-ordinates if the cell is not selected
locn Landscape::getSelectedDistnCell(int dist, int ix) {
	return distns[dist]->getSelectedCell(ix);
}

// Get the dimensions of a specified initial distribution
locn Landscape::getDistnDimensions(int dist) {
	return distns[dist]->getDimensions();
}

// Reset the distribution for a given species so that all cells are deselected
void Landscape::resetDistribution(Species* pSp) {
	// CURRENTLY WORKS FOR FIRST SPECIES ONLY ...
	distns[0]->resetDistribution();
}

//---------------------------------------------------------------------------

// Initialisation cell functions

int Landscape::initCellCount(void) {
	return static_cast<int>(initcells.size());
}

void Landscape::addInitCell(int x, int y) {
	initcells.push_back(new DistCell(x, y));
}

locn Landscape::getInitCell(int ix) {
	return initcells[ix]->getLocn();
}

void Landscape::clearInitCells(void) {
	int ncells = static_cast<int>(initcells.size());
	for (int i = 0; i < ncells; i++) {
		delete initcells[i];
	}
	initcells.clear();
}

//---------------------------------------------------------------------------

// Read landscape file(s)
// Returns error code or zero if read correctly

int Landscape::readLandscape(int fileNum, string habfile, string pchfile, string costfile)
{
	// fileNum == 0 for (first) habitat file and optional patch file
	// fileNum > 0  for subsequent habitat files under the %cover option

#if RS_RCPP
	wstring header;
#else
	string header;
#endif
	int habCode, seq, patchCode, noDataHabCode;
	int noDataPatch = 0;
	int ncols, nrows;
	float habFloat, patchFloat;
	Patch* pPatch;
	simParams sim = paramsSim->getSim();

	if (fileNum < 0) return 19;

#if RS_RCPP
	wifstream ifsHabMap; // habitat file input stream
	wifstream ifsPatchMap; // patch file input stream
#else
	ifstream ifsHabMap; // habitat file input stream
	ifstream ifsPatchMap; // patch file input stream
#endif
	initParams init = paramsInit->getInit();

	// Open habitat file and optionally also patch file
#if RS_RCPP
	ifsHabMap.open(habfile, std::ios::binary);
	if (landraster.utf) {
		// apply BOM-sensitive UTF-16 facet
		ifsHabMap.imbue(std::locale(ifsHabMap.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
	}
#else
	ifsHabMap.open(habfile.c_str());
#endif

	if (!ifsHabMap.is_open()) return 11;

	if (fileNum == 0) {
		if (usesPatches) {
#if RS_RCPP
			ifsPatchMap.open(pchfile, std::ios::binary);
			if (patchraster.utf) {
				// apply BOM-sensitive UTF-16 facet
				ifsPatchMap.imbue(std::locale(ifsPatchMap.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
			}
#else
			ifsPatchMap.open(pchfile.c_str());
#endif
			if (!ifsPatchMap.is_open()) {
				ifsHabMap.close();
				ifsHabMap.clear();
				return 12;
			}
		}
	}

	// read landscape data from header records of habitat file
	// NB headers of all files have already been compared
	ifsHabMap >> header >> ncols 
		>> header >> nrows 
		>> header >> minEast 
		>> header >> minNorth
		>> header >> resol 
		>> header >> noDataHabCode;

#if RS_RCPP
	if (!ifsHabMap.good()) {
		// corrupt file stream
		StreamErrorR(habfile);
		ifsHabMap.close();
		ifsHabMap.clear();
		if (usesPatches) {
			ifsPatchMap.close();
			ifsPatchMap.clear();
		}
		return 131;
	}
#endif

	dimX = ncols; 
	dimY = nrows; 
	minX = maxY = 0; 
	maxX = dimX - 1;
	maxY = dimY - 1;

	if (fileNum == 0) { // First map layer

		// Set initialisation limits to landscape limits
		init.minSeedX = init.minSeedY = 0;
		init.maxSeedX = maxX; 
		init.maxSeedY = maxY;
		paramsInit->setInit(init);

		if (usesPatches) {
			for (int i = 0; i < 5; i++) 
				ifsPatchMap >> header >> patchFloat;
			ifsPatchMap >> header >> noDataPatch;
		}

#if RS_RCPP
		if (!ifsPatchMap.good()) {
			// corrupt file stream
			StreamErrorR(pchfile);
			ifsHabMap.close();
			ifsHabMap.clear();
			ifsPatchMap.close();
			ifsPatchMap.clear();
			return 135;
		}
#endif
		setCellArray();
	}

	// set up bad float values to ensure that valid values are read
	float badHabFloat = -9.0; 
	if (noDataHabCode == -9) badHabFloat = -99.0;
	float badPatchFloat = -9.0; 
	if (noDataPatch == -9) badPatchFloat = -99.0;

	species_id sp = gSingleSpeciesID;
	// Create the matrix patch (even if there is no matrix)
	if (fileNum == 0) {
		patchesList.at(sp).push_back(new Patch(0, 0, sp));
	}
	seq = 1; 			// initial sequential patch landscape
	patchCode = 1; 		// initial patch number for cell-based landscape

	switch (rasterType) {

	case 0: // Raster with habitat codes - 100% habitat each cell
		if (fileNum > 0) return 19; // error condition - should not occur
		
		for (int y = dimY - 1; y >= 0; y--) {
			for (int x = 0; x < dimX; x++) {

				habFloat = badHabFloat;
#if RS_RCPP
				if (ifsHabMap >> habFloat) {
					habCode = (int)habFloat;
					if (usesPatches) {
						patchFloat = badPatchFloat;
						if (ifsPatchMap >> patchFloat) {
							patchCode = (int)patchFloat;
						}
						else { // corrupt file stream
#if !R_CMD
							Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
							StreamErrorR(pchfile);
							ifsHabMap.close();
							ifsHabMap.clear();
							ifsPatchMap.close();
							ifsPatchMap.clear();
							return 132;
						}
					}
				} else { // corrupt file stream
#if !R_CMD
					Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
					StreamErrorR(habfile);
					ifsHabMap.close();
					ifsHabMap.clear();
					if (usesPatches) {
						ifsPatchMap.close();
						ifsPatchMap.clear();
					}
					return 135;
				}
#else
				// Read habitat code in this cell
				ifsHabMap >> habFloat;
				habCode = static_cast<int>(habFloat);
				if (usesPatches) {
					// Read patch code in this cell
					patchFloat = badPatchFloat;
					ifsPatchMap >> patchFloat;
					patchCode = static_cast<int>(patchFloat);
				}
#endif
				if (habCode == noDataHabCode) {
					// No-data cell
					addNewCellToLand(x, y, -1); // x, y is a null cell
				}
				else if (habCode < 0 || (sim.batchMode && (habCode < 1 || habCode > nHabMax))) {
					// Invalid habitat code
#if RS_RCPP && !R_CMD
					Rcpp::Rcout << "Found invalid habitat code." << std::endl;
#endif
					ifsHabMap.close();
					ifsHabMap.clear();
					if (usesPatches) {
						ifsPatchMap.close();
						ifsPatchMap.clear();
					}
					return 13;
				}
				else { // Valid habitat code
					addHabCode(habCode);
					if (usesPatches) {
						if (patchCode < 0 || patchCode == noDataPatch) { // invalid patch code
#if RS_RCPP && !R_CMD
							if (patchCode == noDataPatch) Rcpp::Rcout << "Found patch NA in valid habitat cell." << std::endl;
							else Rcpp::Rcout << "Found negative patch ID in valid habitat cell." << std::endl;
#endif
							ifsHabMap.close();
							ifsHabMap.clear();
							ifsPatchMap.close();
							ifsPatchMap.clear();
							return 14;
						}
						// Does the patch already exists?
						pPatch = patchCode == 0 ? nullptr : ( // matrix cell
							existsPatch(sp, patchCode) ?
							findPatch(sp, patchCode) :
							addNewPatch(sp, seq++, patchCode) // doesn't exist yet: create it
							);
						addNewCellToPatch(pPatch, x, y, habCode);
					}
					else { // cell-based model
						// add cell to landscape (patches created later)
						addNewCellToLand(x, y, habCode);
					}
				}

			} // for x
		} // for y

#if RS_RCPP
		ifsHabMap >> habFloat;
		if (!ifsHabMap.eof()) EOFerrorR(habfile);
		if (usesPatches) {
			ifsPatchMap >> patchFloat;
			if (!ifsPatchMap.eof()) EOFerrorR(pchfile);
		}
#endif
		break;

	case 1: // multiple % cover
		for (int y = dimY - 1; y >= 0; y--) {
			for (int x = 0; x < dimX; x++) {

				habFloat = badHabFloat;
#if RS_RCPP
				if (ifsHabMap >> habFloat) {
					habCode = static_cast<int>(habFloat);
					if (fileNum == 0) { // first habitat cover layer
						if (usesPatches) {
							patchFloat = badPatchFloat;
							if (ifsPatchMap >> patchFloat) {
								patchCode = static_cast<int>(patchFloat);
							}
							else {
								// corrupt file stream
#if !R_CMD
								Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
								StreamErrorR(pchfile);
								ifsHabMap.close();
								ifsHabMap.clear();
								ifsPatchMap.close();
								ifsPatchMap.clear();
								return 135;
							}
						}
#else
				ifsHabMap >> habFloat;
				habCode = static_cast<int>(habFloat);

				if (fileNum == 0) { // first habitat cover layer
					if (usesPatches) {
						patchFloat = badPatchFloat;
						ifsPatchMap >> patchFloat;
						patchCode = static_cast<int>(patchFloat);
						}
#endif
					if (habCode == noDataHabCode) {
						addNewCellToLand(x, y, -1); // add cell only to landscape
					}
					else if (habFloat < 0.0 || habFloat > 100.0) { // invalid cover score
#if RS_RCPP && !R_CMD
						Rcpp::Rcout << "Found invalid habitat cover score." << std::endl;
#endif
						ifsHabMap.close(); 
						ifsHabMap.clear();
						if (usesPatches) {
							ifsPatchMap.close(); 
							ifsPatchMap.clear();
						}
						return 17;
					}
					else {
						if (usesPatches) {
							if (patchCode < 0 || patchCode == noDataPatch) { // invalid patch code
#if RS_RCPP && !R_CMD
								if (patchCode == noDataPatch) Rcpp::Rcout << "Found patch NA in valid habitat cell." << std::endl;
								else Rcpp::Rcout << "Found negative patch ID in valid habitat cell." << std::endl;
#endif
								ifsHabMap.close(); 
								ifsHabMap.clear();
								ifsPatchMap.close();
								ifsPatchMap.clear();
								return 14;
							}
							
							pPatch = patchCode == 0 ? nullptr : ( // matrix cell
								existsPatch(sp, patchCode) ?
								findPatch(sp, patchCode) :
								addNewPatch(sp, seq++, patchCode)
								);
							addNewCellToPatch(pPatch, x, y, habFloat);
						}
						else { // cell-based model
							// add cell to landscape (patches created later)
							addNewCellToLand(x, y, habFloat);
						}
					}
				}
				else { // additional habitat cover layers
					if (habCode != noDataHabCode) {
						if (habFloat < 0.0 || habFloat > 100.0) { // invalid cover score
#if RS_RCPP && !R_CMD
			Rcpp::Rcout << "Found invalid habitat cover score." << std::endl;
#endif
							ifsHabMap.close(); 
							ifsHabMap.clear();
							if (usesPatches) {
								ifsPatchMap.close();
								ifsPatchMap.clear();
							}
							return 17;
						}
						else {
							cells[y][x]->addHabitat(habFloat);
						}
					}
				}
#if RS_RCPP
				else { // not ifsHabMap >> habFloat
					// corrupt file stream
#if RS_RCPP && !R_CMD
					Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
					StreamErrorR(habfile);
					ifsHabMap.close();
					ifsHabMap.clear();
					if (usesPatches) {
						ifsPatchMap.close();
						ifsPatchMap.clear();
					}
					return 133;
				}
#endif

			} // for x
		} // for y

		habsAreIndexed = true; // habitats are already numbered 1...n in correct order

#if RS_RCPP
	ifsHabMap >> habFloat;
	if (!ifsHabMap.eof()) EOFerrorR(habfile);
	if (usesPatches) {
		ifsPatchMap >> patchFloat;
		if (!ifsPatchMap.eof()) EOFerrorR(pchfile);
	}
#endif
		break;

	case 2: // habitat quality
		if (fileNum > 0) return 19; // error condition - should not occur
		for (int y = dimY - 1; y >= 0; y--) {
			for (int x = 0; x < dimX; x++) {

				habFloat = badHabFloat;
#if RS_RCPP
			if (ifsHabMap >> habFloat) {
				habCode = static_cast<int>(habFloat);

			}
			else {
				// corrupt file stream
#if !R_CMD
				Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
				StreamErrorR(habfile);
				ifsHabMap.close();
				ifsHabMap.clear();
				if (usesPatches) {
					ifsPatchMap.close();
					ifsPatchMap.clear();
				}
				return 134;
			}
			if (usesPatches) {
				patchFloat = badPatchFloat;
				if (ifsPatchMap >> patchFloat) {
					patchCode = static_cast<int>(patchFloat);
				}
				else {
					// corrupt file stream
#if RS_RCPP && !R_CMD
					Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
					StreamErrorR(pchfile);
					ifsHabMap.close();
					ifsHabMap.clear();
					ifsPatchMap.close();
					ifsPatchMap.clear();
					return 135;
				}
			}
#else
				ifsHabMap >> habFloat;
				habCode = static_cast<int>(habFloat);
				if (usesPatches) {
					patchFloat = badPatchFloat;
					ifsPatchMap >> patchFloat;
					patchCode = static_cast<int>(patchFloat);
				}
#endif
				if (habCode == noDataHabCode) {
					addNewCellToLand(x, y, -1); // add cell only to landscape
				}
				else if (habFloat < 0.0 || habFloat > 100.0) { // invalid quality score
#if RS_RCPP && !R_CMD
					Rcpp::Rcout << "Found invalid habitat quality score." << std::endl;
#endif
					ifsHabMap.close();
					ifsHabMap.clear();
					if (usesPatches) {
						ifsPatchMap.close();
						ifsPatchMap.clear();
					}
					return 17;
				}
				else {
					if (usesPatches) {
						if (patchCode < 0 || patchCode == noDataPatch) { // invalid patch code
#if RS_RCPP && !R_CMD
							if (patchCode == noDataPatch) Rcpp::Rcout << "Found patch NA in valid habitat cell." << std::endl;
							else Rcpp::Rcout << "Found negative patch ID in valid habitat cell." << std::endl;
#endif
							ifsHabMap.close();
							ifsHabMap.clear();
							ifsPatchMap.close();
							ifsPatchMap.clear();
							return 14;
						}

						if (patchCode != 0) {
							cout << "";
						}
						
						pPatch = patchCode == 0 ? nullptr : ( // matrix cell
							existsPatch(sp, patchCode) ?
							findPatch(sp, patchCode) :
							addNewPatch(sp, seq++, patchCode)
							);
						addNewCellToPatch(pPatch, x, y, habFloat);
					}
					else { // cell-based model
						// add cell to landscape (patches created later)
						addNewCellToLand(x, y, habFloat);
					}
				}

			} // for x
		} // for y

#if RS_RCPP
	ifsHabMap >> habFloat;
	if (!ifsHabMap.eof()) EOFerrorR(habfile);
	if (usesPatches)
	{
		ifsPatchMap >> patchFloat;
		if (!ifsPatchMap.eof()) EOFerrorR(pchfile);
	}
#endif
		break;

	default:
		break;
	} // end switch(rasterType)

	if (ifsHabMap.is_open()) { 
		ifsHabMap.close(); 
		ifsHabMap.clear(); 
	}
	if (ifsPatchMap.is_open()) { 
		ifsPatchMap.close(); 
		ifsPatchMap.clear();
	}

	if (sim.batchMode) {
		if (costfile != "NULL") {
			int retcode = readCosts(costfile);
			if (retcode < 0) return 54;
		}
	}

	return 0;
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

int Landscape::readCosts(string fname)
{

#if RS_RCPP
	wifstream usesCosts; // cost map file input stream
#else
	ifstream costs; // cost map file input stream
#endif

	//int hc,maxYcost,maxXcost,NODATACost,hab;
	int hc, maxYcost, maxXcost, NODATACost;
	float minLongCost, minLatCost; int resolCost;
	float fcost;
#if RS_RCPP
	wstring header;
#else
	string header;
#endif
	Cell* pCell;
#if !RS_RCPP
	simView v = paramsSim->getViews();
#endif

	int maxcost = 0;

	// open cost file
#if RS_RCPP
	usesCosts.open(fname, std::ios::binary);
	if (costsraster.utf) {
		// apply BOM-sensitive UTF-16 facet
		usesCosts.imbue(std::locale(usesCosts.getloc(), new std::codecvt_utf16<wchar_t, 0x10ffff, std::consume_header>));
	}
#else
	costs.open(fname.c_str());
#endif
	// read headers and check that they correspond to the landscape ones
	costs >> header;
#if RS_RCPP
	if (!usesCosts.good()) {
		// corrupt file stream
		StreamErrorR(fname);
		usesCosts.close();
		usesCosts.clear();
		return -181;
	}
	if (header != L"ncols" && header != L"NCOLS") {
#else
	if (header != "ncols" && header != "NCOLS") {
#endif
	costs.close(); costs.clear();
	return -1;
}
double tmpresolCost;
costs >> maxXcost >> header >> maxYcost >> header >> minLongCost;
costs >> header >> minLatCost >> header >> tmpresolCost >> header >> NODATACost;
resolCost = (int) tmpresolCost;

	for (int y = maxYcost - 1; y > -1; y--) {
		for (int x = 0; x < maxXcost; x++) {
#if RS_RCPP
			if (usesCosts >> fcost) {
#else
			costs >> fcost;
#endif
			hc = (int)fcost; // read as float and convert to int
#if RS_RCPP
			}
	else {
		// corrupt file stream
#if RS_RCPP && !R_CMD
		Rcpp::Rcout << "At (x,y) = " << x << "," << y << " :" << std::endl;
#endif
		StreamErrorR(fname);
		usesCosts.close();
		usesCosts.clear();
		return -181;
	}
#endif
			if (hc < 1 && hc != NODATACost) {

#if RS_RCPP && !R_CMD
		Rcpp::Rcout << "Cost map may only contain values of 1 or higher, but found " << fcost << "." << endl;
#endif
		// error - zero / negative cost not allowed
		costs.close(); 
		costs.clear();
		return -999;
	}

	pCell = findCell(x, y);
	if (pCell != nullptr) { // not no-data cell
	    if (hc > 0){ // only if cost value is  above 0 in a data cell
	        pCell->setCost(hc);
		    if (hc > maxcost) maxcost = hc;
	    } else { // if cost value is below 0

#if RS_RCPP && !R_CMD
	    Rcpp::Rcout << "Cost map may only contain values of 1 or higher in habitat cells, but found " << hc << " in cell x: " << x << " y: " << y << "." << endl;
#endif
		throw runtime_error("Found negative- or zero-cost habitat cell.");
	    }

	} // end not no data cell
		}
	}

#if RS_RCPP
	usesCosts >> fcost;
	if (usesCosts.eof()) {
#if !R_CMD
		Rcpp::Rcout << "Costs map loaded." << endl;
#endif
	}
	else EOFerrorR(fname);
#endif

	costs.close(); 
	costs.clear();

	return maxcost;

}

//---------------------------------------------------------------------------

rasterdata CheckRasterFile(string fname)
{
	rasterdata r;
	string header;
	int inint;
	ifstream infile;

	r.ok = true;
	r.errors = r.ncols = r.nrows = r.cellsize = 0;
	r.xllcorner = r.yllcorner = 0.0;

	infile.open(fname.c_str());
	if (infile.is_open()) {
		infile >> header >> r.ncols;
		if (header != "ncols" && header != "NCOLS") r.errors++;
		infile >> header >> r.nrows;
		if (header != "nrows" && header != "NROWS") r.errors++;
		infile >> header >> r.xllcorner;
		if (header != "xllcorner" && header != "XLLCORNER") r.errors++;
		infile >> header >> r.yllcorner;
		if (header != "yllcorner" && header != "YLLCORNER") r.errors++;
		infile >> header >> r.cellsize;
		if (header != "cellsize" && header != "CELLSIZE") r.errors++;
		infile >> header >> inint;
		if (header != "NODATA_value" && header != "NODATA_VALUE") r.errors++;
		infile.close();
		infile.clear();
		if (r.errors > 0) r.ok = false;
	}
	else {
		r.ok = false; 
		r.errors = -111;
	}
	infile.clear();

	return r;
}

//---------------------------------------------------------------------------

// Patch connectivity functions

// Create & initialise connectivity matrix
void Landscape::createConnectMatrix()
{
	for (auto& [sp, connectMatrix] : connectMatrices) {

		if (connectMatrix != nullptr)
			deleteConnectMatrix(sp);

		int npatches = static_cast<int>(patchesList.at(sp).size());
		connectMatrix = new int* [npatches];
		for (int i = 0; i < npatches; i++) {
			connectMatrix[i] = new int[npatches];
			for (int j = 0; j < npatches; j++) 
				connectMatrix[i][j] = 0;
		}
	}
}

// Re-initialise connectivity matrix
void Landscape::resetConnectMatrix() {
	for (auto& [speciesID, connectMatrix] : connectMatrices) {
		if (connectMatrix != nullptr) {
			int npatches = static_cast<int>(patchesList.at(speciesID).size());
			for (int i = 0; i < npatches; i++) {
				for (int j = 0; j < npatches; j++) connectMatrix[i][j] = 0;
			}
		}
	}
}

// Increment connectivity count between two specified patches
void Landscape::incrConnectMatrix(const species_id& speciesID, int originPatchNb, int settlePatchNb) {
	
	int npatches = static_cast<int>(patchesList.at(speciesID).size());
	if (connectMatrices.at(speciesID) == nullptr
		|| originPatchNb < 0 || originPatchNb >= npatches
		|| settlePatchNb < 0 || settlePatchNb >= npatches
		) {
		return;
	}
	connectMatrices.at(speciesID)[originPatchNb][settlePatchNb]++;
}

// Delete connectivity matrix
void Landscape::deleteConnectMatrix(const species_id& speciesID) {
	if (connectMatrices.at(speciesID) != nullptr) {
		int npatches = static_cast<int>(patchesList.at(speciesID).size());
		for (int j = 0; j < npatches; j++) {
			if (connectMatrices.at(speciesID)[j] != nullptr)
				delete connectMatrices.at(speciesID)[j];
		}
		delete[] connectMatrices.at(speciesID);
		connectMatrices.at(speciesID) = nullptr;
	}
}

bool Landscape::closeConnectOfs() {
	if (outConnMat.is_open()) outConnMat.close();
	outConnMat.clear();
	return true;
}

// Write connectivity file headers
bool Landscape::outConnectHeaders()
{
	simParams sim = paramsSim->getSim();

	string name = paramsSim->getDir(2);
	if (sim.batchMode) {
		name += "Batch" + to_string(sim.batchNum) + "_";
		name += "Sim" + to_string(sim.simulation) + "_Land" + to_string(landNum);
	}
	else
		name += "Sim" + to_string(sim.simulation);
	name += "_Connect.txt";
	outConnMat.open(name.c_str());

	outConnMat << "Rep\tYear\tSpecies\tStartPatch\tEndPatch\tNinds" << endl;

	return outConnMat.is_open();
}

#if RS_RCPP
// Write movement paths file headers
void Landscape::outPathsHeaders(int rep, int option)
{
	if (option == -999) { // close the file
		if (outMovePaths.is_open()) outMovePaths.close();
		outMovePaths.clear();
	}
	if (option == 0) { // open the file and write header

		simParams sim = paramsSim->getSim();
		string name = paramsSim->getDir(2);
		if (sim.batchMode) {
			name += "Batch" + to_string(sim.batchNum)
				+ "_Sim" + to_string(sim.simulation)
				+ "_Land" + to_string(landNum)
				+ "_Rep" + to_string(rep);
		}
		else {
			name += "Sim" + to_string(sim.simulation)
				+ "_Rep" + to_string(rep);
		}
		name += "_MovePaths.txt";

		outMovePaths.open(name.c_str());
		if (outMovePaths.is_open()) {
			outMovePaths << "Year\tIndID\tStep\tx\ty\tStatus" << endl;
		}
		else {
			outMovePaths.clear();
		}
	}
}
#endif

void Landscape::outConnect(int rep, int yr) {
	
	int patchnum0, patchnum1;

	for (auto& [speciesID, connectMatrix] : connectMatrices) {

		int npatches = static_cast<int>(patchesList.at(speciesID).size());
		int* emigrants = new int[npatches]; // 1D array to hold emigrants from each patch
		int* immigrants = new int[npatches]; // 1D array to hold immigrants to  each patch

		for (int i = 0; i < npatches; i++) {
			emigrants[i] = immigrants[i] = 0;
		}

		for (int i = 0; i < npatches; i++) {
			patchnum0 = patchesList.at(speciesID)[i]->getPatchNum();
			if (patchnum0 != 0) {
				for (int j = 0; j < npatches; j++) {
					patchnum1 = patchesList.at(speciesID)[j]->getPatchNum();
					if (patchnum1 != 0) {
						emigrants[i] += connectMatrix[i][j];
						immigrants[j] += connectMatrix[i][j];
						if (connectMatrix[i][j] > 0) {
							outConnMat << rep << "\t" << yr
								<< "\t" << patchnum0 << "\t" << patchnum1
								<< "\t" << connectMatrix[i][j] << endl;
						}
					}
				}
			}
		}

		for (int i = 0; i < npatches; i++) {
			patchnum0 = patchesList.at(speciesID)[i]->getPatchNum();
			if (patchnum0 != 0) {
				if (patchesList.at(speciesID)[i]->isSuitable()) {
					outConnMat << rep << "\t" << yr
						<< "\t" << patchnum0 << "\t-999\t" << emigrants[i] << endl;
					outConnMat << rep << "\t" << yr
						<< "\t-999\t" << patchnum0 << "\t" << immigrants[i] << endl;
				}
			}
		}
		delete[] emigrants;
		delete[] immigrants;
	}

}

//---------------------------------------------------------------------------

void Landscape::resetVisits() {
	for (int y = dimY - 1; y >= 0; y--) {
		for (int x = 0; x < dimX; x++) {
			if (cells[y][x] != nullptr) { // not a no-data cell
				cells[y][x]->resetVisits();
			}
		}
	}
}

// Save SMS path visits map to raster text file
void Landscape::outVisits(int rep, int landNr) {

	ofstream outvisits;
	string name;
	simParams sim = paramsSim->getSim();

	if (sim.batchMode) {
		name = paramsSim->getDir(3)
#if RS_RCPP
			+ "Batch" + to_string(sim.batchNum) + "_"
			+ "Sim" + to_string(sim.simulation)
			+ "_Land" + to_string(landNr) + "_Rep" + to_string(rep)
#else
			+ "Batch" + to_string(sim.batchNum) + "_"
			+ "Sim" + to_string(sim.simulation)
			+ "_land" + to_string(landNr) + "_rep" + to_string(rep)
#endif
			+ "_Visits.txt";
	}
	else {
		name = paramsSim->getDir(3)
			+ "Sim" + to_string(sim.simulation)
			+ "_land" + to_string(landNr) + "_rep" + to_string(rep)
			+ "_Visits.txt";
	}
	outvisits.open(name.c_str());

	outvisits << "ncols " << dimX << endl;
	outvisits << "nrows " << dimY << endl;
	outvisits << "xllcorner " << minEast << endl;
	outvisits << "yllcorner " << minNorth << endl;
	outvisits << "cellsize " << resol << endl;
	outvisits << "NODATA_value -9" << endl;

	for (int y = dimY - 1; y >= 0; y--) {
		for (int x = 0; x < dimX; x++) {
			if (cells[y][x] == nullptr) { // no-data cell
				outvisits << "-9 ";
			}
			else {
				outvisits << cells[y][x]->getVisits() << " ";
			}
		}
		outvisits << endl;
	}
	outvisits.close(); 
	outvisits.clear();
}

//---------------------------------------------------------------------------

#ifndef NDEBUG
// Debug only: shortcut setup utilities

landParams createDefaultLandParams(const int& dim) {

	landParams ls_params;
	ls_params.dimX = ls_params.dimY = dim;
	ls_params.minX = ls_params.minY = 0;
	ls_params.maxX = ls_params.maxY = ls_params.dimX - 1;
	ls_params.resol = ls_params.spResol = 1;
	ls_params.rasterType = 0; // habitat types
 
	ls_params.usesPatches = false;
	ls_params.useSpDist = false;
	ls_params.generated = false;
	ls_params.dynamic = false;
	ls_params.landNum = 0;
	ls_params.nHab = ls_params.nHabMax = 0; // irrelevant for habitat codes 
	return ls_params;
}

void testLandscape() {
	// test coordinate system...
}
#endif // NDEBUG

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
