/*
 ============================================================================

 .______    _______     ___   .___________.    __  .___________.
 |   _  \  |   ____|   /   \  |           |   |  | |           |
 |  |_)  | |  |__     /  ^  \ `---|  |----`   |  | `---|  |----`
 |   _  <  |   __|   /  /_\  \    |  |        |  |     |  |     
 |  |_)  | |  |____ /  _____  \   |  |        |  |     |  |     
 |______/  |_______/__/     \__\  |__|        |__|     |__|     
 
 BeatIt - code for cardiovascular simulations
 Copyright (C) 2016 Simone Rossi

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ============================================================================
 */

/*
 * LinearMaterial.cpp
 *
 *  Created on: Oct 21, 2016
 *      Author: srossi
 */

#include "Elasticity/Materials/LinearMaterial.hpp"
#include "libmesh/getpot.h"

namespace BeatIt {

Material* createLinearMaterial()
{
	return new LinearMaterial;
}


LinearMaterial::LinearMaterial()
{
	/*
	 *  The parameters will be stored using the following order
	 *  0: density = rho
	 *  1: shear modulus: mu
	 *  2: bulk modulus: kappa
	 */
	M_parameters.resize(3);

}
LinearMaterial::~LinearMaterial()
{
	// TODO Auto-generated destructor stub
}


void
LinearMaterial::setup(GetPot& data, std::string section)
{
	std::cout << "* LINEAR MATERIAL: Setup. Reading parameters from: " << section << std::endl;
	M_parameters[0] = data(section+"/rho", 0.0); // rho
	double E = data(section+"/E", 0.0);
	double nu = data(section+"/nu", 0.0);
	M_isIncompressible = data(section+"/incompressible", false);
//	M_isIncompressible = false;
	M_parameters[1] = E / ( 2.0 * ( 1 + nu ) );// mu

	if(!M_isIncompressible) M_parameters[2] = E / ( 3.0 * ( 1 - 2 *  nu ) );// kappa
	else M_parameters[2] = 1.0; // trick in order to use the same code as for the compressible case
    std::cout << "\t density = " << M_parameters[0] << std::endl;
    std::cout << "\t shear modulus = " << M_parameters[1] << std::endl;
    std::cout << "\t bulk modulus = " << M_parameters[2] << std::endl;

}

void
LinearMaterial::evaluateStress(ElasticSolverType solverType)
{
    evaluateDeviatoricStress();
    M_total_stress = M_deviatoric_stress;
    switch(solverType)
    {
        case ElasticSolverType::Mixed:
        {
            M_total_stress += M_pressure * M_identity;
            break;
        }
        default:
        case ElasticSolverType::Primal:
        {
            evaluateVolumetricStress();
            M_total_stress += M_volumetric_stress;
            break;
        }
    }
}

void
LinearMaterial::evaluateDeviatoricStress()
{
    double mu = M_parameters[1];
    M_strain =  0.5 * ( M_gradU + M_gradU.transpose() );
    double trEk = M_strain.tr();
//    M_deviatoric_stress = 2.0 * mu * M_strain;
//    if(!M_isIncompressible)     M_deviatoric_stress -= 2.0 * mu / 3.0 * trEk * M_identity;
    M_deviatoric_stress = 2.0 * mu * ( M_strain - trEk / 3.0 * M_identity);

}

void
LinearMaterial::evaluateVolumetricStress()
{
    double trEk = M_strain.tr();
    M_volumetric_stress = M_parameters[2] * trEk * M_identity;
}

void
LinearMaterial::evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q)
{
	double mu = M_parameters[1];
	double kappa = M_parameters[2];
	auto	  dE = 0.5 * (dU + dU.transpose() );
	double trdE = dE.tr();
	M_total_jacobian = 2.0 * mu * ( dE - trdE / 3.0  * M_identity ) + kappa * trdE * M_identity;

}

double
LinearMaterial::evaluatePressure()
{
    double kappa = M_parameters[2];
    M_strain =  0.5 * ( M_gradU + M_gradU.transpose() );
    return kappa * M_strain.tr();
}

double
LinearMaterial::evaluatePressureResidual()
{
//	double res = M_isIncompressible ? 0.0 : M_pressure;
//	std::cout << "res:" <<  res << std::endl;
//    return res -   M_parameters[2] * M_strain.tr();
//    return M_pressure -   M_parameters[2] * M_strain.tr();
    return -   M_parameters[2] * M_strain.tr();
}


void
LinearMaterial::evaluateVolumetricJacobian( const libMesh::TensorValue <double>& dU, double q)
{
	M_volumetric_jacobian = q * M_identity;
}
void
LinearMaterial::evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q)
{
	double mu = M_parameters[1];
	double kappa = M_parameters[2];
	auto	  dE = 0.5 * (dU + dU.transpose() );
	double trdE = dE.tr();
//	M_deviatoric_jacobian = 2.0 * mu * dE;
//	if(!M_isIncompressible) M_deviatoric_jacobian -= 2.0 * mu / 3.0 *  trdE *  M_identity;
		M_deviatoric_jacobian = 2.0 * mu * ( dE - trdE / 3.0 * M_identity);

}

double
LinearMaterial::dpdF(const libMesh::TensorValue <double>&  dF)
{
    return -  M_parameters[2] *  dF.tr();
//    return -  dF.tr();
}



} /* namespace BeatIt */
