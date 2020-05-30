#pragma once

//#include "ComputationEngine.hpp"
//#include "Config.hpp"
//#include "Options.hpp"
//#include "Utilities.hpp"
//#include "ComputationWrapper.hpp"
//#include "FileDescription.hpp"
//#include "InferenceEngine.hpp"
//#include "ParameterManager.hpp"
//#include "OptimizationWrapper.hpp"
//#include "SStruct.hpp"
//#include "Defaults.ipp"
#include <string>
class Options;

void InitOptions(Options &options);

float pfunc(char* sequence, char* constraints = "?");
char* predict_struct(char* sequence, char* constraints = "?");
float c_fold(char* sequence, char* structure);
float c_energy_of_structure(char* sequence, char* structure);
