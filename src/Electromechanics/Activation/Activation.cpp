/*
 * Activation.cpp
 *
 *  Created on: Jul 13, 2017
 *      Author: srossi
 */

#include "Electromechanics/Activation/Activation.hpp"

namespace BeatIt
{

Activation::Activation(int numVar, const std::string name)
: M_numVariables(numVar)
, M_variablesNames(numVar)
, M_modelName(name)
{
}

Activation::~Activation()
{
}




} /* namespace BeatIt */
