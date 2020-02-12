/*
 * IonicModel.cpp
 *
 *  Created on: Aug 6, 2016
 *      Author: Simone Rossi
 */

#include "Electrophysiology/IonicModels/IonicModel.hpp"
#include "libmesh/getpot.h"


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
  , M_membrane_capacitance(1.0)
{
}


void IonicModel::setup(GetPot& data, std::string section)
{
    M_membrane_capacitance = data(section+"/Cm", 1.0 ); //uF/cm^2
    M_surface_to_volume_ratio = data(section+"/Chi", 1400.0 );// 1/cm
}

double
IonicModel::evaluateIonicCurrentH(std::vector<double>& variables, double appliedCurrent, double dt, double h)
{
	return evaluateIonicCurrent(variables, appliedCurrent, dt);
}


void
IonicModel::updateVariables(std::vector<double>& variables, std::vector<double>& /*rhs*/, double appliedCurrent, double dt, bool /*overwrite*/)
{
    updateVariables(variables, appliedCurrent, dt);
};


void
IonicModel::solve( std::vector<double>& variables,
                   double appliedCurrent,
                   double dt)
{
    updateVariables(variables, appliedCurrent, dt);
    // evaluateIonicCurrent does not containe appliedCurrent
    // Cm dV/dt = - Iion - Istim
    variables[0] += dt * (- evaluateIonicCurrent(variables, appliedCurrent, dt) - appliedCurrent);
}



} // namespace BeatIt
