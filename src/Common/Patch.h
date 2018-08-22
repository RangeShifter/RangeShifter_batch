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
Bocedi G., Palmer S.C.F., Pe�er G., Heikkinen R.K., Matsinos Y.G., Watts K.
and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
eco-evolutionary dynamics and species� responses to environmental changes.
Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

Authors: Greta Bocedi & Steve Palmer, University of Aberdeen

Last updated: 17 March 2017 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef PatchH
#define PatchH

#include <vector>
using namespace std;

#include "Version.h"
#include "Parameters.h"
#include "Cell.h"
#include "Species.h"

#if VCL
#include <VCLTee.Chart.hpp>
#endif

//---------------------------------------------------------------------------

struct patchLimits {
	int xMin,xMax,yMin,yMax;
};
struct patchPopn {
	intptr pSp,pPop; // pointers to Species and Population cast as integers
};

class Patch{
public:
	Patch(
		int,	 // internal sequential number
		int		 // patch id number
	);
	~Patch();
	int getSeqNum(void);
	int getPatchNum(void);
	int getNCells(void);
	patchLimits getLimits(void); // Returns the minimum and maximum co-ordinates of the patch
	bool withinLimits( // Does the patch fall (partially) within a specified rectangle?
		patchLimits // structure holding the SW and NE co-ordinates of the rectangle
	);
	void resetLimits(void); // Reset minimum and maximum co-ordinates of the patch
	void addCell(
		Cell*,	// pointer to the Cell to be added to the Patch
		int,int	// x (column) and y (row) co-ordinates of the Cell
	);
	locn getCell( // Return co-ordinates of a specified cell
		int			// index no. of the Cell within the vector cells
	);
	locn getCentroid(void); // Return co-ordinates of patch centroid
	void removeCell(
		Cell*	// pointer to the Cell to be removed from the Patch
	);
	Cell* getRandomCell(void);
	void setSubComm(
		intptr		// pointer to the Sub-community cast as an integer
	);
	intptr getSubComm(void);
	void addPopn(
		patchPopn // structure holding pointers to Species and Population cast as integers
	);
	intptr getPopn( // return pointer (cast as integer) to the Population of the Species
		intptr // pointer to Species cast as integer
	);
	void resetPopn(void);
	void resetPossSettlers(void);
	void incrPossSettler( // Record the presence of a potential settler within the Patch
		Species*,	// pointer to the Species
		int				// sex of the settler
	);
	int getPossSettlers( // Get number of a potential settlers within the Patch
		Species*, // pointer to the Species
		int       // sex of the settlers
	);
	void setCarryingCapacity( // Calculate total Patch carrying capacity (no. of inds)
		Species*, 		// pointer to the Species
		patchLimits,	// current min and max limits of landscape
		float,				// global stochasticity value (epsilon) for the current year
		short,				// no. of habitat classes in the Landscape
		short,				// rasterType (see Landscape)
		short,				// landscape change index (always zero if not dynamic)
		bool					// TRUE if there is a gradient in carrying capacity across the Landscape
	);
	float getK(void);
#if VCL
	// for GUI version, draw the Patch on the screen
	void drawCells(TCanvas*,float,int,rgb);
#else
	// dummy function for batch version
	void drawCells(float,int,rgb);
#endif

	private:
	int patchSeqNum;// sequential patch number - patch 0 is reserved for the inter-patch matrix
	int patchNum; 	// patch number as supplied by the user (not forced to be sequential)
	int nCells;			// no. of cells in the patch
	int xMin,xMax,yMin,yMax; 	// min and max cell co-ordinates
	int x,y;				// centroid co-ordinates (approx.)
	intptr subCommPtr; // pointer (cast as integer) to sub-community associated with the patch
	// NOTE: FOR MULTI-SPECIES MODEL, PATCH WILL NEED TO STORE K FOR EACH SPECIES
	float localK;		// patch carrying capacity (individuals)
	bool changed;
// NOTE: THE FOLLOWING ARRAY WILL NEED TO BE MADE SPECIES-SPECIFIC...
	short nTemp[NSEXES];						// no. of potential settlers in each sex

	std::vector <Cell*> cells;
	std::vector <patchPopn> popns;

};

//---------------------------------------------------------------------------

extern paramStoch *paramsStoch;
#if RSRANDOM
extern RSrandom *pRandom;
#else
extern StochasticLib1 *pRandom;
#endif

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

#endif