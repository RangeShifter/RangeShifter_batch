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

 RangeShifter v2.0 SubCommunity

 Implements the SubCommunity class

 There is ONE instance of a SubCommunity for each Patch in the Landscape
 (including the matrix). The SubCommunity holds a number of Populations, one for
 each Species represented in the simulation.
 CURRENTLY the number of Populations withn a SubCommunity is LIMITED TO ONE.

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

 Last updated: 25 June 2021 by Greta Bocedi

 ------------------------------------------------------------------------------*/

#ifndef SubCommunityH
#define SubCommunityH

#include <vector>
#include <algorithm>
using namespace std;

#include "Parameters.h"
#include "Landscape.h"
#include "Population.h"

//---------------------------------------------------------------------------

class SubCommunity {

public:
	SubCommunity(Patch* pPch, int subCommId);
	~SubCommunity();
private:
	int subCommNum;	// SubCommunity number
		// 0 is reserved for the SubCommunity in the inter-patch matrix
	Patch *pPatch;
	std::vector <Population*> popns;
};

extern paramGrad* paramsGrad;
extern paramStoch* paramsStoch;
extern paramInit* paramsInit;

//---------------------------------------------------------------------------
#endif
