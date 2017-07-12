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

#include "Elasticity/Materials/IsotropicMaterial.hpp"
#include "libmesh/getpot.h"
namespace BeatIt {

Material* createIsotropicMaterial()
{
	return new IsotropicMaterial;
}

IsotropicMaterial::IsotropicMaterial()
{
	/*
	 *  The parameters will be stored using the following order
	 *  0: density = rho
	 *  1: shear modulus: mu
	 *  2: bulk modulus: kappa
	 */
	M_parameters.resize(3);
	M_FA = M_identity;
    M_FAinv = M_identity;
    M_CAinv = M_identity;

    // For debugging purposes
    // M_FA = M_identity;
    // M_FA(0,0) -= 0.2;
    // M_FA(1,1) += 0.25;
    // M_FAinv = M_FA.inverse();
    // M_CAinv = M_FAinv * M_FAinv;

}

IsotropicMaterial::~IsotropicMaterial() {
	// TODO Auto-generated destructor stub
}

void
IsotropicMaterial::setup(GetPot& data, std::string section)
{
	std::cout << "* ISOTROPIC MATERIAL: Setup. Reading parameters from: " << section << std::endl;
	M_parameters[0] = data(section+"/rho", 0.0); // rho

	std::string type = data(section+"/type", "UNKNOWN");
	std::map<std::string, IsotropicCL> cl_map;
	cl_map["neohookean"] = IsotropicCL::Neohookean;
    cl_map["mooneyrivlin"] = IsotropicCL::MooneyRivlin;
    cl_map["exponential"] = IsotropicCL::Exponential;
    cl_map["fung"] = IsotropicCL::Fung;
    auto it_cl = cl_map.find(type);
    if( it_cl != cl_map.end() ) M_cl = it_cl->second;
    else
    {
        std::cout << "--- IsotropicMaterial: type " << type << " is not known" << std::endl;
        throw std::runtime_error("IsotropicMaterial setup error");
    }

    M_active = data(section+"/active", false);

    M_isIncompressible = true;
    switch(M_cl)
    {
        case IsotropicCL::Neohookean:
        {
            double E = data(section+"/E", 0.0);
            double nu = data(section+"/nu", 0.0);
            M_parameters[1] = E / ( 2.0 * ( 1 + nu ) );// mu
            M_parameters[2] = 1.0;// kappa
            M_kappa = M_parameters[2];
            std::cout << "\t density = " << M_parameters[0] << std::endl;
            std::cout << "\t shear modulus = " << M_parameters[1] << std::endl;
            std::cout << "\t bulk modulus = " << M_parameters[2] << std::endl;
            M_tau = 0.5 /  M_parameters[1];

            break;
        }
        case IsotropicCL::MooneyRivlin:
        {
            double c1 = data(section+"/c1", 0.0);
            double c2 = data(section+"/c2", 0.0);
            M_parameters[1] = c1;// c1
            M_parameters[2] = c2;// c2
            M_kappa = data(section+"/k", 1.0);
            std::cout << "\t density = " << M_parameters[0] << std::endl;
            std::cout << "\t c1 = " << M_parameters[1] << std::endl;
            std::cout << "\t c2 = " << M_parameters[2] << std::endl;
            M_tau = 0.5 / (c1 + c2 );

            break;
        }
        case IsotropicCL::Exponential:
        {
            double a = data(section+"/a", 0.0);
            double b = data(section+"/b", 0.0);
            M_parameters[1] = a;// a
            M_parameters[2] = b;// v
            M_kappa = data(section+"/k", 1.0);

            std::cout << "\t density = " << M_parameters[0] << std::endl;
            std::cout << "\t a = " << M_parameters[1] << std::endl;
            std::cout << "\t b = " << M_parameters[2] << std::endl;
            M_tau = 0.5 /  a;

            break;
        }
        case IsotropicCL::Fung:
        {
            double C = data(section+"/C", 0.0);
            double b = data(section+"/b", 0.0);
            M_parameters[1] = C;// C
            M_parameters[2] = b;// b
            M_kappa = data(section+"/k", 1.0);

            std::cout << "\t density = " << M_parameters[0] << std::endl;
            std::cout << "\t C = " << M_parameters[1] << std::endl;
            std::cout << "\t b = " << M_parameters[2] << std::endl;
            M_tau = 0.5 /  C;

            break;
        }
        default:
        {
            break;
        }
    }
    M_density = M_parameters[0];
    double tau_coeff = data(section+"/tau", 1.0); // stabilization coefficient
    M_tau *= tau_coeff; // rho

}


void
IsotropicMaterial::evaluateDeviatoricStress()
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
    M_W1  = W1 (M_I1, M_I2, M_Jk);
    M_W11 = W11(M_I1, M_I2, M_Jk);
    M_W12 = W12(M_I1, M_I2, M_Jk);
    M_W2  = W2 (M_I1, M_I2, M_Jk);
    M_W21 = W21(M_I1, M_I2, M_Jk);
    M_W22 = W22(M_I1, M_I2, M_Jk);

    if(M_active)
    {
        M_deviatoric_stress  = 2.0 * ( M_W1 + M_W2 * M_I1 ) * ( M_CAinv - M_I1 / 3.0 * M_Cinvk);
        M_deviatoric_stress -= 2.0 * M_W2 * ( M_Ck * M_CAinv * M_CAinv - M_J2 / 3.0 * M_Cinvk);
    }
    else
    {
        M_deviatoric_stress  = 2.0 * ( M_W1 + M_W2 * M_I1 ) * ( M_identity - M_I1 / 3.0 * M_Cinvk);
        M_deviatoric_stress -= 2.0 * M_W2 * ( M_Ck - M_J2 / 3.0 * M_Cinvk);
    }
}

void
IsotropicMaterial::evaluateVolumetricStress()
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
	double kappa = M_kappa;

	double p = kappa * ( M_Jk - 1);
	M_volumetric_stress = M_Jk * p * M_Cinvk;
}

void
IsotropicMaterial::updateVariables()
{
    M_Fk = M_identity + M_gradU;
    MaterialUtilities::cof(M_Fk, M_Hk);
    M_Ck = M_Fk.transpose() * M_Fk;
    M_Cinvk = M_Ck.inverse();
    M_Jk = M_Fk.det();
    M_I3 = M_Jk * M_Jk;

    if(M_active)
    {
        M_I1 = M_Ck.contract(M_CAinv);
        auto tmp = M_Ck * M_CAinv;
        M_J2 = tmp.contract(tmp);
        M_I2 = 0.5 * (M_I1 * M_I1 - M_J2 );
    }
    else
    {
        M_I1 = M_Ck.tr();
        M_J2 = M_Ck.contract(M_Ck);
        M_I2 = 0.5 * (M_I1 * M_I1 - M_J2 );
    }
}

void
IsotropicMaterial::evaluateStress(ElasticSolverType solverType)
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

	evaluateDeviatoricStress();
    M_total_stress = M_deviatoric_stress;
    switch(solverType)
    {
        case ElasticSolverType::Mixed:
        {
            M_total_stress += M_pressure * M_Cinvk;
            break;
        }
        default:
        case ElasticSolverType::Primal:
        {
            std::cout << "--- ISOTROPIC MATERIAL: Primal Formulation not coded for isotropic material!!!" << std::endl;
            throw std::runtime_error("Primal Formulation not coded for isotropic material");
            break;
        }
    }
	M_PK1 = M_Fk * M_total_stress;

//	double mu = M_parameters[1];
//	M_PK1 = mu * (M_Fk - M_identity) +  M_Jk * M_pressure * M_Cinvk;
}


void
IsotropicMaterial::evaluateVolumetricJacobian( const libMesh::TensorValue <double>& dU, double q)
{
	M_volumetric_jacobian = q * M_Fk * M_Cinvk;
}


void
IsotropicMaterial::evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q )
{
	auto  dF = dU;
	auto  dC = M_Fk.transpose() * dF + dF.transpose() * M_Fk;
    auto CinvdCCinv = M_Cinvk * dC * M_Cinvk;
    double coeff_1 = M_W11 + M_W12 * M_I1 + M_W2 + M_I1 * M_W21 + M_I1 * M_I1 * M_W22;
    double coeff_2 = M_W12 + M_I1 * M_W22;
    auto Jac = 0.0 * M_identity;
    if(M_active)
    {
        double CAinvdC = M_CAinv.contract(dC);
        auto CCainv2 = M_Ck * M_CAinv * M_CAinv;
        double CCainv2dC = CCainv2.contract(dC);

        double IdC = M_identity.contract(dC);
        double CdC = M_Ck.contract(dC);
        auto dW1 = ( M_W11 + M_W12 * M_I1 ) * CAinvdC - M_W12 * CCainv2dC;
        auto dW2 = ( M_W21 + M_W22 * M_I1 ) * CAinvdC - M_W22 * CCainv2dC;
        auto dI1 = CAinvdC;
        auto dJ2 = 2.0 * CCainv2dC;
        // delta S1
        Jac  = 2.0 * ( dW1 + dW2 * M_I1 + M_W2 * dI1) * ( M_CAinv - M_I1/ 3.0 * M_Cinvk);
        Jac += 2.0 * ( M_W1 + M_W2 * M_I1 ) * ( CinvdCCinv * M_I1 / 3.0 - dI1 / 3.0 * M_Cinvk);
        // delta S1
        Jac -= 2.0 * dW2  * ( CCainv2 - M_J2 / 3.0 * M_Cinvk);
        Jac -= 2.0 * M_W2 * ( dC * M_CAinv * M_CAinv + CinvdCCinv * M_J2 / 3.0 - dJ2 / 3.0 * M_Cinvk);
    }
    else
    {

        double IdC = M_identity.contract(dC);
        double CdC = M_Ck.contract(dC);
        auto dW1 = ( M_W11 + M_W12 * M_I1 ) * IdC - M_W12 * CdC;
        auto dW2 = ( M_W21 + M_W22 * M_I1 ) * IdC - M_W22 * CdC;
        auto dI1 = IdC;
        auto dJ2 = 2.0 * CdC;
        // delta S1
        Jac  = 2.0 * ( dW1 + dW2 * M_I1 + M_W2 * dI1) * ( M_identity - M_I1/ 3.0 * M_Cinvk);
        Jac += 2.0 * ( M_W1 + M_W2 * M_I1 ) * ( CinvdCCinv * M_I1 / 3.0 - dI1 / 3.0 * M_Cinvk);
        // delta S1
        Jac -= 2.0 * dW2  * ( M_Ck - M_J2 / 3.0 * M_Cinvk);
        Jac -= 2.0 * M_W2 * ( dC + CinvdCCinv * M_J2 / 3.0 - dJ2 / 3.0 * M_Cinvk);
    }
	/*
    *  dS4 = ( J p )' * dJ/dC * dC
    *             -   J p C^-1 dC C^-1
    *         =  0.5 * J * ( C^-1 : dC ) ( p  ) * C^-1
    *             -   J p C^-1 dC C^-1
    */
    double p = M_pressure;
	Jac -= p * CinvdCCinv;
	M_deviatoric_jacobian = dU * M_total_stress + M_Fk * Jac;
}

void
IsotropicMaterial::evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q)
{
    std::cout << "--- ISOTROPIC MATERIAL: evaluateJacobian not coded for isotropic material!!!" << std::endl;
    throw std::runtime_error("Primal Formulation not coded for isotropic material");
}


double
IsotropicMaterial::evaluatePressure()
{
    double kappa = M_parameters[2];
    M_Fk = M_identity + M_gradU;
    M_Jk = M_Fk.det();
    return kappa * ( M_Jk - 1);
}

double
IsotropicMaterial::evaluatePressureResidual()
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
	    std::cout << "--- ISOTROPIC MATERIAL: evaluatePressureResidual compressible case not coded for isotropic material!!!" << std::endl;
	    throw std::runtime_error("Primal Formulation not coded for isotropic material");
	}
	return - RHS;

}

double
IsotropicMaterial::dpdF(const libMesh::TensorValue <double>&  dF)
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
        std::cout << "--- ISOTROPIC MATERIAL: dpdF compressible case not coded for isotropic material!!!" << std::endl;
        throw std::runtime_error("Primal Formulation not coded for isotropic material");
		return 0.0;
	}
}

double
IsotropicMaterial::d2U( double J)
{
    // For the time being this material is only incompressible
    // we return 1 to be compatible with the code
    return 1.0;
}
double
IsotropicMaterial::d3U( double J)
{
    return 0.0;
}


/////// MATERIALS:
/////// CONSTITUTIVE LAWS:



double
IsotropicMaterial::W1 (double I1, double I2, double J)
{
    double W1 = 0.0;
    switch(M_cl)
    {
        case IsotropicCL::Neohookean:
        case IsotropicCL::MooneyRivlin:
        {
            // mu = M_parameters[1]
            W1 = 0.5 * M_parameters[1];
            break;
        }
        case IsotropicCL::Exponential:
        {
            // a = M_parameters[1]
            // b = M_parameters[2]
            W1 = 0.5 * M_parameters[1] * std::exp(M_parameters[2] * (I1 - 3) );
            break;
        }
        case IsotropicCL::Fung:
        {
            // C = M_parameters[1]
            // b = M_parameters[2]
        	double Q = 0.25 * M_parameters[2] * ( I1 * I1 - 2.0 * I2 - 2.0 * I1 + 3 );
        	double dQdI1 = 0.25 * M_parameters[2] * ( 2.0 * I1 - 2.0 );
        	W1 = 0.5 * M_parameters[1] * dQdI1 * std::exp(Q);
            break;
        }
        default:
        {
            break;
        }
    }
    return W1;
}
double
IsotropicMaterial::W11(double I1, double I2, double J)
{
    double W11 = 0.0;
    switch(M_cl)
    {
        case IsotropicCL::Exponential:
        {
            // a = M_parameters[1]
            // b = M_parameters[1]
            W11 = 0.5 * M_parameters[1] * M_parameters[2] * std::exp(M_parameters[2] * (I1 - 3) );
            break;
        }
        case IsotropicCL::Fung:
        {
            // C = M_parameters[1]
            // b = M_parameters[2]
        	double Q = 0.25 * M_parameters[2] * ( I1 * I1 - 2.0 * I2 - 2.0 * I1 + 3 );
        	double dQdI1 = 0.25 * M_parameters[2] * ( 2.0 * I1 - 2.0 );
        	double d2QdI1 = 0.25 * M_parameters[2] * ( 2.0 );
        	W11 = 0.5 * M_parameters[1] * ( dQdI1 * dQdI1 + d2QdI1 ) * std::exp(Q);
            break;
        }

        case IsotropicCL::Neohookean:
        case IsotropicCL::MooneyRivlin:
        default:
        {
            break;
        }
    }
    return W11;
}

double
IsotropicMaterial::W12(double I1, double I2, double J)
{
    double W12 = 0.0;
    switch(M_cl)
    {
        case IsotropicCL::Fung:
        {
            // C = M_parameters[1]
            // b = M_parameters[2]
        	double Q = 0.25 * M_parameters[2] * ( I1 * I1 - 2.0 * I2 - 2.0 * I1 + 3 );
        	double dQdI1 = 0.25 * M_parameters[2] * ( 2.0 * I1 - 2.0 );
        	double dQdI2 = 0.25 * M_parameters[2] * ( - 2.0 );
        	//double d2QdI1dI2 = 0.0;
        	//W12 = M_parameters[1] * ( dQdI1 * dQdI2 + d2QdI1dI2 ) * std::exp(Q);
        	W12 = 0.5 * M_parameters[1] * ( dQdI1 * dQdI2 ) * std::exp(Q);
            break;
        }
        case IsotropicCL::Exponential:
        case IsotropicCL::Neohookean:
        case IsotropicCL::MooneyRivlin:
        default:
        {
            break;
        }
    }
    return W12;
}

double
IsotropicMaterial::W2 (double I1, double I2, double J)
{
    double W2 = 0.0;
    switch(M_cl)
    {
        case IsotropicCL::MooneyRivlin:
        {
            // mu = M_parameters[1]
            W2 = 0.5 * M_parameters[1];
            break;
        }
        case IsotropicCL::Fung:
        {
            // C = M_parameters[1]
            // b = M_parameters[2]
        	double Q = 0.25 * M_parameters[2] * ( I1 * I1 - 2.0 * I2 - 2.0 * I1 + 3 );
        	double dQdI2 = 0.25 * M_parameters[2] * ( - 2.0 );
        	W2 = 0.5 * M_parameters[1] * dQdI2 * std::exp(Q);
            break;
        }
        case IsotropicCL::Neohookean:
        case IsotropicCL::Exponential:
        default:
        {
            break;
        }
    }
    return W2;
}

double
IsotropicMaterial::W21(double I1, double I2, double J)
{
    double W21 = 0.0;
    switch(M_cl)
    {
        case IsotropicCL::Fung:
        {
            // C = M_parameters[1]
            // b = M_parameters[2]
        	double Q = 0.25 * M_parameters[2] * ( I1 * I1 - 2.0 * I2 - 2.0 * I1 + 3 );
        	double dQdI1 = 0.25 * M_parameters[2] * ( 2.0 * I1 - 2.0 );
        	double dQdI2 = 0.25 * M_parameters[2] * ( - 2.0 );
        	//double d2QdI1dI2 = 0.0;
        	//W21 = M_parameters[1] * ( dQdI1 * dQdI2 + d2QdI1dI2 ) * std::exp(Q);
        	W21 = 0.5 * M_parameters[1] * ( dQdI1 * dQdI2 ) * std::exp(Q);
            break;
        }
        case IsotropicCL::Exponential:
        case IsotropicCL::Neohookean:
        case IsotropicCL::MooneyRivlin:
        default:
        {
            break;
        }
    }
    return W21;
}

double
IsotropicMaterial::W22(double I1, double I2, double J)
{
    double W22 = 0.0;
    switch(M_cl)
    {
        case IsotropicCL::Fung:
        {
            // C = M_parameters[1]
            // b = M_parameters[2]
        	double Q = 0.25 * M_parameters[2] * ( I1 * I1 - 2.0 * I2 - 2.0 * I1 + 3 );
        	double dQdI2 = 0.25 * M_parameters[2] * ( - 2.0 );
        	//double d2QdI2 = 0.0;
        	//W22 = M_parameters[1] * ( dQdI2 * dQdI2 + d2QdI2 ) * std::exp(Q);
        	W22 = 0.5 * M_parameters[1] * ( dQdI2 * dQdI2 ) * std::exp(Q);
            break;
        }
        case IsotropicCL::Exponential:
        case IsotropicCL::Neohookean:
        case IsotropicCL::MooneyRivlin:
        default:
        {
            break;
        }
    }
    return W22;
}


} /* namespace BeatIt */



