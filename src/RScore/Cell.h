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

 RangeShifter v2.0 Cell

 Implements the following classes:

 Cell - Landscape cell

 DistCell - Initial species distribution cell

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 14 January 2021 by Steve Palmer

 ------------------------------------------------------------------------------*/

#ifndef CellH
#define CellH

#include <vector>
using namespace std;

#include "Parameters.h"
#include "Species.h"

#ifdef _OPENMP
#include <atomic>
#include <mutex>
#endif

class Patch;

//---------------------------------------------------------------------------

class Patch; // Forward-declaration of the Patch class

struct array3x3f { float cell[3][3]; }; 	// neighbourhood cell array (SMS)
struct smscosts { int cost; array3x3f* effcosts; };	// cell costs for SMS

// Landscape cell

class Cell {
public:
	// Constructor for habitat codes
	Cell(int xx, int yy, int hab, set<species_id> spLabels);
	// Constructor for habitat % cover or habitat quality
	Cell(int xx, int yy, float hab, set<species_id> spLabels);

	~Cell();
	void addHabIndex(short hx);
	void changeHabIndex(short ix, short hx);
	int getHabIndex(int ix);
	int nHabitats();
	void addHabitat(float q); // habitat prop or cell quality score
	float getHabitat(int ix); // Get habitat proportion / quality score
	void setPatch(species_id whichSpecies, Patch* p);
	Patch* getPatch(species_id whichSpecies);
	locn getLocn();
	void setEnvDev(float devVal);
	float getEnvDev();
	void updateEps(float ac, float randpart); // Update local environmental stochasticity (epsilon)
	float getEps();
	void setCost(species_id sp, int costVal);
	int getCost(species_id sp);
	void resetCost();
	array3x3f getEffCosts(species_id sp);
	void setEffCosts(species_id sp, array3x3f effCostsNeighbours);
	void resetEffCosts(species_id sp); // Reset the effective cost, but not the cost, of the cell
	void resetVisits();
	void incrVisits(species_id sp);
	unsigned long int getVisits(species_id sp);
#ifdef _OPENMP
	std::unique_lock<std::mutex> lockCost(void);
#endif
	void declareOverlappingPatches() const;
private:
	int x, y;		// cell co-ordinates

	map<species_id, Patch*> patches; // which patch the cell belongs to for each species

	float envDev;	// local environmental deviation (static, in range -1.0 to +1.0)
	float eps;		// local environmental stochasticity (epsilon) (dynamic, from N(0,std))
#ifdef _OPENMP
	map<species_id, std::atomic<unsigned long int>> visits; // no. of times the cell is visited by each species
#else
	map<species_id, unsigned long int> visits; // no. of times the cell is visited by each species
#endif
	map<species_id, smscosts*> smsData;

	vector <short> habIxx; 	// habitat indices (rasterType=0)
							// initially, habitat codes are loaded, then converted to index nos.
							// once landscape is fully loaded
	vector <float> habitats;	// habitat proportions (rasterType=1) or quality (rasterType=2)

#ifdef _OPENMP
	std::mutex cost_mutex;
#endif
};

//---------------------------------------------------------------------------

// Initial species distribution cell

class DistCell {
public:
	DistCell(
		int,	// x co-ordinate
		int		// y co-ordinate
	);
	~DistCell();
	void setCell(
		bool	// TRUE if cell is to be initialised, FALSE if not
	);
	bool toInitialise(
		locn	// structure holding co-ordinates of cell
	);
	bool selected(void);
	locn getLocn(void);

private:
	int x, y;					// cell co-ordinates
	bool initialise;  // cell is to be initialised
};

//---------------------------------------------------------------------------

#endif
