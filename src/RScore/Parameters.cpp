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

#include "Parameters.h"

// Environmental stochasticity parameters

paramStoch::paramStoch() {
	stoch = false; 
	local = false; 
	inK = false; 
	usesLocalExt = false;
	ac = 0.0; std = 0.25;
	locExtProb = 0.1;
}

paramStoch::~paramStoch() {}

void paramStoch::setStoch(envStochParams e)
{
	stoch = e.usesStoch; 
	local = e.stochIsLocal; 
	inK = e.inK; 
	usesLocalExt = e.usesLocalExt;
	if (e.ac >= 0.0 && e.ac < 1.0) ac = e.ac;
	if (e.std > 0.0 && e.std <= 1.0) std = e.std;
	locExtProb = e.locExtProb;
}

bool paramStoch::envStoch() { return stoch; }

envStochParams paramStoch::getStoch()
{
	envStochParams e;
	e.usesStoch = stoch; 
	e.stochIsLocal = local; 
	e.inK = inK; 
	e.usesLocalExt = usesLocalExt;
	e.ac = ac; 
	e.std = std;
	e.locExtProb = locExtProb;
	return e;
}

//---------------------------------------------------------------------------

// Simulation parameters

paramSim::paramSim(const string& pathToProjDir) :
	dir{pathToProjDir}
{
	simulation = 0;
	reps = years = 1;
	batchMode = absorbing = false;
}

paramSim::~paramSim() { }

void paramSim::setSim(simParams s) {
	if (s.batchNum >= 0) batchNum = s.batchNum;
	if (s.simulation >= 0) simulation = s.simulation;
	if (s.reps >= 1) reps = s.reps;
	if (s.years >= 1) years = s.years;
	batchMode = s.batchMode;
	absorbing = s.absorbing;
	fixReplicateSeed = s.fixReplicateSeed;
}

void paramSim::setGeneticSim(string patchSamplingOption, bool outputGeneticValues, bool outputWeirCockerham, bool outputWeirHill, int outputStartGenetics, int outputGeneticInterval) {
	this->patchSamplingOption = patchSamplingOption;
	this->outputGenes = outputGeneticValues;
	this->outputWeirCockerham = outputWeirCockerham;
	this->outputWeirHill = outputWeirHill;
	this->outputStartGenetics = outputStartGenetics;
	this->outputGeneticInterval = outputGeneticInterval;
}

simParams paramSim::getSim() {
	simParams s;
	s.batchNum = batchNum;
	s.simulation = simulation; s.reps = reps; s.years = years;
	s.outRange = outRange; s.outOccup = outOccup; s.outPop = outPop; s.outInds = outInds;
	s.outTraitsCells = outTraitsCells; s.outTraitsRows = outTraitsRows; s.outConnect = outConnect;
	s.outStartPop = outStartPop; s.outStartInd = outStartInd;
	s.outStartTraitCell = outStartTraitCell; s.outStartTraitRow = outStartTraitRow;
	s.outStartConn = outStartConn;
	s.outIntRange = outIntRange;
	s.outIntOcc = outIntOcc; s.outIntPop = outIntPop;
	s.outIntInd = outIntInd;
	s.outIntTraitCell = outIntTraitCell;
	s.outIntTraitRow = outIntTraitRow;
	s.outIntConn = outIntConn;
	s.batchMode = batchMode;
	s.absorbing = absorbing;
	s.traitInt = traitInt;
#if RS_RCPP
	s.outStartPaths = outStartPaths;
	s.outIntPaths = outIntPaths;
	s.outPaths = outPaths;
	s.ReturnPopRaster = ReturnPopRaster;
	s.CreatePopFile = CreatePopFile;
#endif
	s.patchSamplingOption = patchSamplingOption;
	s.outputGeneValues = outputGenes;
	s.outputWeirCockerham = outputWeirCockerham;
	s.outputWeirHill = outputWeirHill;
	s.outStartGenetics = outputStartGenetics;
	s.outputGeneticInterval = outputGeneticInterval;

	return s;
}

int paramSim::getSimNum() { return simulation; }

// return directory name depending on option specified
string paramSim::getDir(int option) {
	string s;
	switch (option) {
	case 0: // working directory
		s = dir;
		break;
#if LINUX_CLUSTER || RS_RCPP
	case 1: // Inputs folder
		s = dir + "Inputs/";
		break;
	case 2: // Outputs folder
		s = dir + "Outputs/";
		break;
	case 3: // Maps folder
		s = dir + "Output_Maps/";
		break;
#else
	case 1: // Inputs folder
		s = dir + "Inputs\\";
		break;
	case 2: // Outputs folder
		s = dir + "Outputs\\";
		break;
	case 3: // Maps folder
		s = dir + "Output_Maps\\";
		break;
#endif
	default:
		s = "ERROR_ERROR_ERROR";
	}
	return s;
}

string to_string(const TraitType& tr) {
	switch (tr)
	{
	case NEUTRAL: return "NEUTRAL";
	case GENETIC_LOAD: return "GENETIC_LOAD";
	case GENETIC_LOAD1: return "GENETIC_LOAD1";
	case GENETIC_LOAD2: return "GENETIC_LOAD2";
	case GENETIC_LOAD3: return "GENETIC_LOAD3";
	case GENETIC_LOAD4: return "GENETIC_LOAD4";
	case GENETIC_LOAD5: return "GENETIC_LOAD5";

	case E_D0: return "E_D0";
	case E_D0_M: return "E_D0_M";
	case E_D0_F: return "E_D0_F";
	case E_ALPHA: return "E_ALPHA";
	case E_ALPHA_M: return "E_ALPHA_M";
	case E_ALPHA_F: return "E_ALPHA_F";
	case E_BETA: return "E_BETA";
	case E_BETA_M: return "E_BETA_M";
	case E_BETA_F: return "E_BETA_F";

	case S_S0: return "S_S0";
	case S_S0_M: return "S_S0_M";
	case S_S0_F: return "S_S0_F";
	case S_ALPHA: return "S_ALPHA";
	case S_ALPHA_M: return "S_ALPHA_M";
	case S_ALPHA_F: return "S_ALPHA_F";
	case S_BETA: return "S_BETA";
	case S_BETA_M: return "S_BETA_M";
	case S_BETA_F: return "S_BETA_F";

	case CRW_STEPLENGTH: return "CRW_STEPLENGTH";
	case CRW_STEPCORRELATION: return "CRW_STEPCORRELATION";
	case KERNEL_MEANDIST_1: return "KERNEL_MEANDIST_1";
	case KERNEL_MEANDIST_2: return "KERNEL_MEANDIST_2";
	case KERNEL_MEANDIST_1_F: return "KERNEL_MEANDIST_1_F";
	case KERNEL_MEANDIST_2_F: return "KERNEL_MEANDIST_2_F";
	case KERNEL_MEANDIST_1_M: return "KERNEL_MEANDIST_1_M";
	case KERNEL_MEANDIST_2_M: return "KERNEL_MEANDIST_2_M";
	case KERNEL_PROBABILITY: return "KERNEL_PROBABILITY";
	case KERNEL_PROBABILITY_F: return "KERNEL_PROBABILITY_F";
	case KERNEL_PROBABILITY_M: return "KERNEL_PROBABILITY_M";

	case SMS_DP: return "SMS_DP";
	case SMS_GB: return "SMS_GB";
	case SMS_ALPHADB: return "SMS_ALPHADB";
	case SMS_BETADB: return "SMS_BETADB";
	case INVALID_TRAIT: return "INVALID_TRAIT";
	default: return "";
	}
}

string to_string(const GenParamType& param) {
	switch (param)
	{
	case MEAN: return "MEAN";
	case SD: return "SD";
	case MIN: return "MIN";
	case MAX: return "MAX";
	case SHAPE: return "SHAPE";
	case SCALE: return "SCALE";
	case INVALID: return "INVALID";
	default: return "";
	}
}

string to_string(const DistributionType& dist) {
	switch (dist)
	{
	case UNIFORM: return "UNIFORM";
	case NORMAL: return "NORMAL";
	case GAMMA: return "GAMMA";
	case NEGEXP: return "NEGEXP";
	case SCALED: return "SCALED";
	case KAM: return "KAM";
	case SSM: return "SSM";
	case NONE: return "NONE";
	default: return "";
	}
}

string to_string(const ExpressionType& expr) {
	switch (expr)
	{
	case AVERAGE: return "AVERAGE";
	case ADDITIVE: return "ADDITIVE";
	case NOTEXPR: return "NOTEXPR";
	case MULTIPLICATIVE: return "MULTIPLICATIVE";
	default: return "";
	}
}

#if RS_RCPP
bool paramSim::getReturnPopRaster(void) { return ReturnPopRaster; }
bool paramSim::getCreatePopFile(void) { return CreatePopFile; }
#endif

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

