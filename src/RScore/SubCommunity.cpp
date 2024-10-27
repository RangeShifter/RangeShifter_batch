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

#include "SubCommunity.h"
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

SubCommunity::SubCommunity(Patch* pPch, int subCommId) {
	subCommNum = subCommId;
	pPatch = pPch;
	// record the new sub-community no. in the patch
	pPatch->setSubComm(this);
}

SubCommunity::~SubCommunity() {
	pPatch->setSubComm(0);
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		delete popns[i];
	}
	popns.clear();
}

//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


