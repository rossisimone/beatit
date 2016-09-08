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
                        int numGatingVar,
                        const std::string& name,
                        CellType cell_type)
  : M_numVariables(numVar)
  , M_numGatingVariables(numGatingVar)
  , M_cellType(cell_type)
  , M_variablesNames(numVar-1)
  , M_ionicModelName(name)
{
}

} // namespace BeatIt
