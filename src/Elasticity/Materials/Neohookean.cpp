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
 * IsotropicMaterial.cpp
 *
 *  Created on: Oct 21, 2016
 *      Author: srossi
 */

#include "Elasticity/Materials/Neohookean.hpp"
#include "libmesh/getpot.h"
namespace BeatIt {

Material* createNeohookean()
{
	return new Neohookean;
}

Neohookean::Neohookean()
{
	/*
	 *  The parameters will be stored using the following order
	 *  0: density = rho
	 *  1: shear modulus: mu
	 *  2: bulk modulus: kappa
	 */
	M_parameters.resize(3);
	M_volumetricEnergy = U::Quadratic;
}

Neohookean::~Neohookean() {
	// TODO Auto-generated destructor stub
}

void
Neohookean::setup(GetPot& data, std::string section)
{
	std::cout << "* ISOTROPIC MATERIAL: Setup. Reading parameters from: " << section << std::endl;
	M_parameters[0] = data(section+"/rho", 0.0); // rho
	double E = data(section+"/E", 0.0);
	double nu = data(section+"/nu", 0.0);
	M_parameters[1] = E / ( 2.0 * ( 1 + nu ) );// mu

	if( nu < 0.5)
    {
	    M_isIncompressible = false;
	    M_parameters[2] = E / ( 3.0 * ( 1 - 2 *  nu ) );// kappa
	    std::string u = data(section+"/U", "quadratic"); // rho
	    if(u == "liu")
	    {
	        M_volumetricEnergy = U::Liu;
	    }
	    else if(u == "reese")
	    {
            M_volumetricEnergy = U::Reese;
	    }

    }
	else
	{
        M_isIncompressible = true;
        // trick
        M_parameters[2] = 1.0;// kappa

	}
    std::cout << "\t density = " << M_parameters[0] << std::endl;
    std::cout << "\t shear modulus = " << M_parameters[1] << std::endl;
    std::cout << "\t bulk modulus = " << M_parameters[2] << std::endl;
    if( M_volumetricEnergy == U::Liu )
    {
        std::cout << "\t using Liu volumetric energy " << std::endl;
    }
    else
    {
        std::cout << "\t using Quadratic volumetric energy  = " << M_parameters[2] << std::endl;
    }
    M_density = M_parameters[0];

    double ctau = data(section+"/ctau",1.0); // rho

    M_tau = ctau * 0.5 /  M_parameters[1];
    std::cout << "\t tau = " << M_tau << std::endl;
}


void
Neohookean::evaluateDeviatoricStress()
{
    /*
     *  Consider the energy W(I1bar, I2bar) + U(J)
     *
     *  TODO: add the I2bar terms
     *
     *  Consider only W(I1bar), with Fbar = J^(-1/3) * F
     *
     *  Then
     *
     *  W1 = dW / dI1bar
     *
     *  Pdev = W1 * dI1bar / dF = 2 W1 J^(-2/3) ( F - I1 / 3 * F^-T )
     *
     *  P = F * S
     *  Sdev = 2 * W1 * J^(-2/3) * ( I - I1 / 3 * C^-1)
     */
	double mu = M_parameters[1];
	double W1 = mu / 2.0;

	double Jm23 = std::pow(M_Jk, -2.0/3.0);
	if(M_volumetricEnergy == U::Reese)
	{
	    Jm23 = 1.0;
	}

	double I1 = M_Ck.tr();
	M_deviatoric_stress = 2.0 * W1 * Jm23 * ( M_identity - I1 / 3.0 * M_Cinvk);
}

void
Neohookean::evaluateVolumetricStress()
{
    /*
     *  Consider the energy W(I1bar, I2bar) + U(J)
     *
     *  U = k * U(J)
     *  Pvol = U'(J) * J * F^-T
     *  p = U'(J)
     *
     *  P = F * S
     *  Svol =  J p C^-1
     */
    double p = evaluatePressure();
	M_volumetric_stress = M_Jk * p * M_Cinvk;
}


void
Neohookean::evaluateStress(ElasticSolverType solverType)
{
    /*
     *  Consider the energy W(I1bar, I2bar) + U(J)
     *
     *  TODO: add the I2bar terms
     *
     *  Consider only W(I1bar) + U(J), with Fbar = J^(-1/3) * F
     *
     *  Then
     *
     *  W1 = dW / dI1bar
     *
     *  Pdev = W1 * dI1bar / dF = 2 W1 J^(-2/3) ( F - I1 / 3 * F^-T )
     *
     *  U = k * U(J)
     *  Pvol = U'(J) * J * F^-T
     *  p = U'(J)
     *
     *  P = F * S
     *  S = 2 * W1 * J^(-2/3) * ( I - I1 / 3 * C^-1) + J p C^-1
     */


	M_Fk = M_identity + M_gradU;
	MaterialUtilities::cof(M_Fk, M_Hk);
	M_Ck = M_Fk.transpose() * M_Fk;
	M_Cinvk = M_Ck.inverse();
	M_Jk = M_Fk.det();


	evaluateDeviatoricStress();
    M_total_stress = M_deviatoric_stress;
    switch(solverType)
    {
        case ElasticSolverType::Mixed:
        {
            M_total_stress += M_Jk * M_pressure * M_Cinvk;
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
	M_PK1 = M_Fk * M_total_stress;

//	double mu = M_parameters[1];
//	M_PK1 = mu * (M_Fk - M_identity) +  M_Jk * M_pressure * M_Cinvk;
}


void
Neohookean::evaluateVolumetricJacobian( const libMesh::TensorValue <double>& dU, double q)
{
	M_volumetric_jacobian = q *  M_Jk * M_Fk * M_Cinvk;
}


void
Neohookean::evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q )
{
    /*
     *  Consider the energy W(I1bar, I2bar) + U(J)
     *
     *  TODO: add the I2bar terms
     *
     *  Consider only W(I1bar) + U(J), with Fbar = J^(-1/3) * F
     *
     *  Then
     *
     *  W1 = dW / dI1bar
     *
     *  p
     *
     *  S = 2 * W1 * J^(-2/3) * ( I - I1 / 3 * C^-1) + J p C^-1
     *
     *  dP = d(FS) = dF * S + F * dS
     *
     *  dS = dS1 + dS2 + dS3 + dS4
     *
     *  W11 = dW1 / dI1bar
     *  dS1 = 2 * ( 2 * W11 * J^(-2/3) ( dC - I1 / 3 C^-1:dC ) ) * J^(-2/3) * ( I - I1 / 3 * C^-1)
     *
     *  dJ/dC  = 0.5 * J * C^-1 : dC
     *  dJ^(-2/3)/dC = -2 / 3 * J^(-5/3) dJ/dC = - 1 / 3 * J^(-2/3) * C^-1 : dC
     *  dS2 = - 2 * W1 * dJ^(-2/3)/dC * ( I - I1 / 3 * C^-1)
     *
     *  dI1 = I : dC
     *  dS3 = - 2 * W1 * J^(-2/3) * ( I:dC ) / 3 * C^-1
     *        + 2 * W1 * J^(-2/3) * I1 / 3 * C^-1 * dC * C^-1
     *
     *  dS4 = ( J p )' * dJ/dC * dC - J p C^-1 dC C^-1
     *      =  0.5 * J * ( C^-1 : dC ) ( p ) * C^-1    -      J p C^-1 dC C^-1
     */
	double mu = M_parameters[1];
	double kappa = M_parameters[2];
	double W1 = mu / 2.0;

	double Jm23 = std::pow(M_Jk, -2.0/3.0);
    if(M_volumetricEnergy == U::Reese)
    {
        Jm23 = 1.0;
    }

	auto  dF = dU;
	auto  dC = M_Fk.transpose() * dF + dF.transpose() * M_Fk;
	double I1 = M_Ck.tr();
	/*
     *  dS1 = 2 * ( 2 * W11 * J^(-2/3) ( dC - I1 / 3 C^-1:dC ) ) * J^(-2/3) * ( I - I1 / 3 * C^-1)
     *
     *  For the neohookean W11 = 0.0
	 */

    /*
     *  dJ/dC  = 0.5 * J * C^-1 : dC
     *  dJ^(-2/3)/dC = -2 / 3 * J^(-5/3) dJ/dC = - 1 / 3 * J^(-2/3) * C^-1 : dC
     *  dS2 = 2 * W1 * dJ^(-2/3)/dC * ( I - I1 / 3 * C^-1)
     */
    auto dJm23dC = - Jm23 / 3.0 * M_Cinvk.contract(dC);
    if(M_volumetricEnergy == U::Reese)
    {
        dJm23dC = 0.0;
    }

	auto Jac = 2.0 * W1 * dJm23dC * ( M_identity - I1 / 3.0 * M_Cinvk);

	/*
     *  dS3 = - 2 * W1 * J^(-2/3) * ( I:dC ) / 3 * C^-1
     *             + 2 * W1 * J^(-2/3) *  I1 / 3 * C^-1 * dC * C^-1
     */
	Jac -= 2.0 * W1 * Jm23 * dC.tr() / 3.0 * M_Cinvk;
	Jac += 2.0 * W1 * Jm23 * I1 / 3.0 * M_Cinvk * dC * M_Cinvk;
	M_deviatoric_jacobian = dU * M_deviatoric_stress + M_Fk * Jac;

	/*
    *  dS4 = ( J p )' * dJ/dC * dC
    *             -   J p C^-1 dC C^-1
    *         =  0.5 * J * ( C^-1 : dC ) ( p  ) * C^-1
    *             -   J p C^-1 dC C^-1
    */
    double p = M_pressure;
	Jac += 0.5 *  p * M_Jk * M_Cinvk.contract(dC) * M_Cinvk;
	Jac -= M_Jk * p  * M_Cinvk * dC * M_Cinvk;
	M_deviatoric_jacobian = dU * M_total_stress + M_Fk * Jac;
}

void
Neohookean::evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q)
{
    /*
     *  Consider the energy W(I1bar, I2bar) + U(J)
     *
     *  TODO: add the I2bar terms
     *
     *  Consider only W(I1bar) + U(J), with Fbar = J^(-1/3) * F
     *
     *  Then
     *
     *  W1 = dW / dI1bar
     *
     *  p = U'(J)
     *
     *  S = 2 * W1 * J^(-2/3) * ( I - I1 / 3 * C^-1) + J p C^-1
     *
     *  dP = d(FS) = dF * S + F * dS
     *
     *  dS = dS1 + dS2 + dS3 + dS4
     *
     *  W11 = dW1 / dI1bar
     *  dS1 = 2 * ( 2 * W11 * J^(-2/3) ( dC - I1 / 3 C^-1:dC ) ) * J^(-2/3) * ( I - I1 / 3 * C^-1)
     *
     *  dJ/dC  = 0.5 * J * C^-1 : dC
     *  dJ^(-2/3)/dC = -2 / 3 * J^(-5/3) dJ/dC = - 1 / 3 * J^(-2/3) * C^-1 : dC
     *  dS2 = - 2 * W1 * dJ^(-2/3)/dC * ( I - I1 / 3 * C^-1)
     *
     *  dI1 = I : dC
     *  dS3 = - 2 * W1 * J^(-2/3) * ( I:dC ) / 3 * C^-1
     *        + 2 * W1 * J^(-2/3) * I1 / 3 * C^-1 * dC * C^-1
     *
     *  dS3 = ( J U'(J) )' * dJ/dC * dC - J p C^-1 dC C^-1
     *      =  0.5 * J * ( C^-1 : dC ) ( U'(J) + J U''(J) ) * C^-1 - J p C^-1 dC C^-1
     */
	double mu = M_parameters[1];
	double kappa = M_parameters[2];
	double W1 = mu / 2.0;

	double Jm23 = std::pow(M_Jk, -2.0/3.0);
    if(M_volumetricEnergy == U::Reese)
    {
        Jm23 = 1.0;
    }
	auto  dF = dU;
	auto  dC = M_Fk.transpose() * dF + dF.transpose() * M_Fk;
	double I1 = M_Ck.tr();
	/*
     *  dS1 = 2 * ( 2 * W11 * J^(-2/3) ( dC - I1 / 3 C^-1:dC ) ) * J^(-2/3) * ( I - I1 / 3 * C^-1)
     *
     *  For the neohookean W11 = 0.0
	 */

    /*
     *  dJ/dC  = 0.5 * J * C^-1 : dC
     *  dJ^(-2/3)/dC = -2 / 3 * J^(-5/3) dJ/dC = - 1 / 3 * J^(-2/3) * C^-1 : dC
     *  dS2 = 2 * W1 * dJ^(-2/3)/dC * ( I - I1 / 3 * C^-1)
     */
    auto dJm23dC = - Jm23 / 3.0 * M_Cinvk.contract(dC);
    if(M_volumetricEnergy == U::Reese)
    {
        dJm23dC = 0.0;
    }
	auto Jac = 2.0 * W1 * dJm23dC * ( M_identity - I1 / 3.0 * M_Cinvk);

	/*
     *  dS3 = - 2 * W1 * J^(-2/3) * ( I:dC ) / 3 * C^-1
     *        + 2 * W1 * J^(-2/3) * I1 / 3 * C^-1 * dC * C^-1
     */
	Jac -= 2.0 * W1 * Jm23 * M_identity.contract(dC) / 3.0 * M_Cinvk;
	Jac += 2.0 * W1 * Jm23 * I1 / 3.0 * M_Cinvk * dC * M_Cinvk;

	/*
    *  dS4 = ( J p )' * dJ/dC * dC - J p C^-1 dC C^-1
    *      =  0.5 * J * ( C^-1 : dC ) ( U'(J) + J U''(J) ) * C^-1 - J p C^-1 dC C^-1
    */
    double p = evaluatePressure();
    double dJp = p + M_Jk * d2U(M_Jk);
	Jac += 0.5 * dJp * M_Jk * M_Cinvk.contract(dC) * M_Cinvk;
	Jac -= M_Jk * p  * M_Cinvk * dC * M_Cinvk;
	M_total_jacobian = dU * M_total_stress + M_Fk * Jac;

}


double
Neohookean::evaluatePressure()
{
    double kappa = M_parameters[2];
    M_Fk = M_identity + M_gradU;
    M_Jk = M_Fk.det();
    double dU = 0.0;
    switch(M_volumetricEnergy)
    {
        case U::Liu:
        {
            dU = std::log(M_Jk);
            break;
        }
        case U::Reese:
        {
            dU = 0.5 * (M_Jk - 1.0 / M_Jk );
            break;
        }
        case U::Quadratic:
        default:
        {
            dU = ( M_Jk - 1);
            break;
        }
    }
    return kappa * dU;
}

double
Neohookean::evaluatePressureResidual()
{
	/*
	 * 	Evaluate the RHS of the pressure equation with a minus, that is
	 * 	p = RHS
	 * 	p - RHS = 0
	 * 	return -RHS
	 */
	double RHS;
	if(M_isIncompressible)
	{
		RHS = ( M_Jk - 1);
	}
	else
	{
		RHS = evaluatePressure();
	}
	return - RHS;

}

double
Neohookean::dpdF(const libMesh::TensorValue <double>&  dF)
{
	/*
	 * 	 p - k U'(J) = 0
	 *
	 * 	 dp / dF = k U''(J) dJ/DF = k U''(J) J F^-T
	 *
	 * 	dPvol/dp =
	 */
	auto cofF = M_Jk * M_Fk * M_Cinvk;
	if(M_isIncompressible)
	{
		return - cofF.contract(dF);
	}
	else
	{
		return - d2U(M_Jk) * cofF.contract(dF);
	}
}

double
Neohookean::d2U( double J)
{
    double kappa = M_parameters[2];
    double d2U = 0.0;
    switch(M_volumetricEnergy)
    {
        case U::Liu:
        {
            d2U = 1.0 / M_Jk;
            break;
        }
        case U::Reese:
        {
            d2U = 0.5 * (1.0 + 1.0 / M_Jk / M_Jk );
            break;
        }

        case U::Quadratic:
        default:
        {
            d2U = 1.0;
            break;
        }
    }
    return kappa * d2U;
}
double
Neohookean::d3U( double J)
{
    double kappa = M_parameters[2];
    double d3U = 0.0;
    switch(M_volumetricEnergy)
    {
        case U::Liu:
        {
            d3U =-1.0 / M_Jk  / M_Jk;
            break;
        }
        case U::Reese:
        {
            d3U = - 0.5 / M_Jk / M_Jk / M_Jk;
            break;
        }
        case U::Quadratic:
        default:
        {
            break;
        }
    }
    return kappa * d3U;
}




} /* namespace BeatIt */



