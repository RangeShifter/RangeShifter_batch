/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2026 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Roslyn Henry, Théo Pannetier, Jette Wolff, Damaris Zurell
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
 * File Created by Jette Wolff
 --------------------------------------------------------------------------
 */


/*------------------------------------------------------------------------------

 RangeShifter v2.0 Parameters

 Implements the following classes:

 paramManagement  - Management parameters
 paramTranslocation  - Translocation parameters


 Last updated: 12 March 2024 by Jette Reeg

 ------------------------------------------------------------------------------
 */

#ifndef ManagementH
#define ManagementH

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <map>
#if RS_RCPP
#include <RcppArmadillo.h>
#endif
#include "Parameters.h"
#include "Species.h"
#include "Cell.h"
#include "Landscape.h"
#include "Population.h"
#include "Community.h"
using namespace std;

//---------------------------------------------------------------------------

/*
 * Management settings
*/

// Structure for management parameters
struct managementParams {
    bool usesTranslocation; // Translocation
};

// Structure for translocation parameters
struct translocationParams {
    double catching_rate; // Catching rate
    std::vector<int> translocation_years; // Number of years of translocation -> will be increased at the beginning of a simulation
    std::map< int, species_id> species; // Which one species is being translocated
    std::map< int, std::vector <locn> > source; // Source patch or cell: should be a vector of arrays
    std::map< int, std::vector <locn> > target; // Target patch or cell
    std::map< int, std::vector <int> > nb; // number of ttanslocated individuals
    std::map< int, std::vector <int> > min_age; // Minimum age of translocated individuals
    std::map< int, std::vector <int> > max_age; // Maximum age of translocated individuals
    std::map< int, std::vector <int> > stage; // Stage of translocated individuals
    std::map< int, std::vector <int> > sex; // Sex of translocated individuals
};

//---------------------------------------------------------------------------

class Management {
public:
    Management();
    ~Management();
    void setManagementParams( // function to set management parameters
       const managementParams	// structure holding general management parameters
    );
    managementParams getManagementParams(); // get management parameters
    void setTranslocationParams( // function to set translocation parameters
            const translocationParams	// structure holding translocation parameters
    );
    translocationParams getTranslocationParams();
    void translocate(int yr, Landscape* pLandscape, const speciesMap_t& allSpecies, Community* pComm);
    bool isTranslocationYear(const int yr) {
        return (this->usesTranslocation &&
            std::find(translocation_years.begin(), translocation_years.end(), yr) != translocation_years.end());
    }

    bool usesTranslocation; // Translocation
    double catching_rate; // Catching rate
    bool non_dispersed; // whether non-dispersed individuals should be translocated
    std::vector<int> translocation_years; // Number of years of translocation -> should be a dynamic vector
    std::map< int, species_id > species; // Which one species is being translocated
    std::map< int, std::vector <locn> > source; // Source patch or cell: should be a vector of arrays
    std::map< int, std::vector <locn> > target; // Target patch or cell
    std::map< int, std::vector <int> > nb; // number of ttanslocated individuals
    std::map< int, std::vector <int> > min_age; // Minimum age of translocated individuals
    std::map< int, std::vector <int> > max_age; // Maximum age of translocated individuals
    std::map< int, std::vector <int> > stage; // Stage of translocated individuals
    std::map< int, std::vector <int> > sex; // Sex of translocated individuals
};

//---------------------------------------------------------------------------

extern paramSim *paramsSim;

//---------------------------------------------------------------------------
#endif // MANAGEMENTH
