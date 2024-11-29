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

class Patch;

//---------------------------------------------------------------------------

struct array3x3f { float cell[3][3]; }; 	// neighbourhood cell array (SMS)
struct smscosts { int cost; array3x3f* effcosts; };	// cell costs for SMS

// Landscape cell

class Cell {
public:
	// Constructor for habitat codes
	Cell(int xx, int yy, Patch* patch, int hab);

	// Constructor for habitat % cover or habitat quality
	Cell(int xx, int yy, Patch* patch, float hab);

	~Cell();
	void addHabIndex(short hx);
	void changeHabIndex(short ix, short hx);
	int getHabIndex(int ix);
	int nHabitats(void);
	void addHabitat(float q); // habitat prop or cell quality score
	float getHabitat(int ix); // Get habitat proportion / quality score
	void setPatch(Patch* p);
	Patch* getPatch();
	locn getLocn();
	void setEnvDev(
		float	// local environmental deviation
	);
	float getEnvDev();
	void setEnvVal(float e);
	float getEnvVal();
	void updateEps(float ac, float randpart); // Update local environmental stochasticity (epsilon)
	float getEps();
	void setCost(
		int		// cost value for SMS
	);
	int getCost();
	void resetCost();
	array3x3f getEffCosts();
	void setEffCosts(
		array3x3f	// 3 x 3 array of effective costs for neighbouring cells
	);
	void resetEffCosts(); // Reset the effective cost, but not the cost, of the cell
	void resetVisits();
	void incrVisits();
	unsigned long int getVisits();

private:
	int x, y;			// cell co-ordinates
	Patch* pPatch; 	// pointer (cast as integer) to the Patch to which cell belongs
	// NOTE: THE FOLLOWING ENVIRONMENTAL VARIABLES COULD BE COMBINED IN A STRUCTURE
	// AND ACCESSED BY A POINTER ...
	float envVal; // environmental value, representing one of:
	// gradient in K, r or extinction probability
	float envDev;	// local environmental deviation (static, in range -1.0 to +1.0)
	float eps;		// local environmental stochasticity (epsilon) (dynamic, from N(0,std))
	unsigned long int visits; // no. of times square is visited by dispersers
	smscosts* smsData;

	vector <short> habIxx; 		// habitat indices (rasterType=0)
	// NB initially, habitat codes are loaded, then converted to index nos.
	//    once landscape is fully loaded
	vector <float> habitats;	// habitat proportions (rasterType=1) or quality (rasterType=2)
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
