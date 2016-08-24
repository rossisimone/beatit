/*
 * IonicModel.cpp
 *
 *  Created on: Aug 6, 2016
 *      Author: Simone Rossi
 */

#include "Electrophysiology/IonicModels/IonicModel.hpp"

namespace BeatIt
{


IonicModel::IonicModel( int numVar,
		                int numGatingVar)
  : M_numVariables(numVar)
  , M_numGatingVariables(numGatingVar)
  , M_cellType(CellType::Endocardial)
  , M_variablesNames(numVar-1)
{
}

} // namespace BeatIt
