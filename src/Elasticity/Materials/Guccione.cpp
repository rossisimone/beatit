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
 * Guccione.cpp
 *
 *  Created on: Oct 21, 2016
 *      Author: srossi
 */

#include "Elasticity/Materials/Guccione.hpp"
#include "libmesh/getpot.h"

namespace BeatIt {

Material* createGuccione()
{
	return new Guccione;
}

Guccione::Guccione()
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

Guccione::~Guccione() {
	// TODO Auto-generated destructor stub
}

void
Guccione::setup(GetPot& data, std::string section)
{
	std::cout << "* GUCCIONE MATERIAL: Setup. Reading parameters from: " << section << std::endl;
	M_parameters.resize(6);
	M_parameters[0] = data(section+"/rho", 0.0); // rho
    M_active = data(section+"/active", false);

    double nu = data(section+"/nu", -1.0);
    if(nu <= -1 || nu >= 0.5) M_isIncompressible = true;
    else M_isIncompressible = false;

    M_C = data(section+"/C", 0.0);
    M_bff = data(section+"/bff", 0.0);
    M_bfs = data(section+"/bfs", 0.0);
    M_bt = data(section+"/bt", 0.0);
    M_parameters[1] = M_C;
    M_parameters[2] = M_bff;
    M_parameters[3] = M_bfs;
    M_parameters[4] = M_bt;
    // double nu = data(section+"/nu", 0.0);
    //M_parameters[1] = E / ( 2.0 * ( 1 + nu ) );// mu
    if(!M_isIncompressible)
        M_kappa = data(section+"/k", 1.0);
    else
        M_kappa=1.0;
    M_parameters[5] = M_kappa;

    std::cout << "\t density = " << M_parameters[0] << std::endl;
    std::cout << "\t C = " << M_C << std::endl;
    std::cout << "\t bff = " << M_bff << std::endl;
    std::cout << "\t bt = " << M_bt << std::endl;
    std::cout << "\t bfs = " << M_bfs << std::endl;
    std::cout << "\t bulk modulus = " << M_kappa << std::endl;
    std::cout << "\t nu = " << nu << std::endl;
    M_tau = 0.5 /  M_parameters[1];

    M_density = M_parameters[0];

    // stabilization
    M_tau = 0.5 / M_C;
    double tau_coeff = data(section+"/ctau", 1.0); // stabilization coefficient
    M_tau *= tau_coeff; // rho

    M_A(0,0) = M_bff;
    //
    M_A(1,1) = M_bt;
    M_A(2,2) = M_bt;
    M_A(1,2) = M_bt;
    M_A(2,1) = M_bt;

    //
    M_A(0,1) = M_bfs;
    M_A(0,2) = M_bfs;
    M_A(1,0) = M_bfs;
    M_A(2,0) = M_bfs;


}


void
Guccione::evaluateDeviatoricStress()
{
    /*

     */
    M_EFk = M_R.transpose() * M_Ek * M_R;
    for(int idim = 0; idim < 3; idim++)
        for(int jdim = 0; jdim < 3; jdim++)
            M_AEFk(idim,jdim) = M_A(idim,jdim)*M_EFk(idim,jdim);

    M_Q = M_AEFk.contract(M_EFk);
    M_expQ = std::exp(M_Q);
    M_CexpQ = M_C * M_expQ;

    // values for loops:
    if(M_active)
    {
        std::cout << "Active Guccione Not Yet Coded" << std::endl;
        throw std::runtime_error("error");
    }
    else
    {
        M_Sk = M_CexpQ * M_R * M_AEFk * M_R.transpose();
        //M_deviatoric_stress = M_Sk - M_Sk.contract(M_Ck) / 3.0 * M_Cinvk;
        M_deviatoric_stress = M_Sk;
    }
}

void
Guccione::evaluateVolumetricStress()
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

	//double p = kappa * std::log(M_Jk);
    double p = kappa * (M_Jk - 1.0);
	M_volumetric_stress = M_Jk * p * M_Cinvk;
}

void
Guccione::updateVariables()
{
    M_Fk = M_identity + M_gradU;
    MaterialUtilities::cof(M_Fk, M_Hk);
    M_Ck = M_Fk.transpose() * M_Fk;
    M_Cinvk = M_Ck.inverse();
    M_Jk = M_Fk.det();

    M_Ek = 0.5 * (M_Ck-M_identity);

    if(M_Jk <= 0)
    {
        std::cout << "Jk is negative: " << M_Jk << std::endl;
        M_Fk.print(std::cout);
        exit(1);
    }
    //std::cout << "Jk: " << M_Jk << std::endl;
    M_I3 = M_Jk * M_Jk;

    if(M_active)
    {
        std::cout << "Active Guccione Not Yet Coded" << std::endl;
        throw std::runtime_error("error");
    }
    else
    {
        M_n0 = M_f0.cross(M_s0);
        for(int idim = 0; idim < 3; idim++)
        {
            M_R(idim, 0) = M_f0(idim);
            M_R(idim, 1) = M_s0(idim);
            M_R(idim, 2) = M_n0(idim);
        }
    }
}

void
Guccione::evaluateStress(ElasticSolverType solverType)
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
    switch(solverType)
    {
        case ElasticSolverType::Mixed:
        {
            M_total_stress = M_deviatoric_stress;
            M_total_stress += M_pressure * M_Cinvk;
            break;
        }
        default:
        case ElasticSolverType::Primal:
        {
            M_total_stress = M_Sk;
            evaluateVolumetricStress();
            M_total_stress += M_volumetric_stress;

            //std::cout << "--- ISOTROPIC MATERIAL: Primal Formulation not coded for isotropic material!!!" << std::endl;
            //throw std::runtime_error("Primal Formulation not coded for isotropic material");
            break;
        }
    }
	M_PK1 = M_Fk * M_total_stress;

//	double mu = M_parameters[1];
//	M_PK1 = mu * (M_Fk - M_identity) +  M_Jk * M_pressure * M_Cinvk;
}


void
Guccione::evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q )
{
	auto  dF = dU;
	auto  dC = M_Fk.transpose() * dF + dF.transpose() * M_Fk;
    auto  dE = 0.5 * dC;
    auto CinvdCCinv = M_Cinvk * dC * M_Cinvk;
    //auto Jac = 0.0 * M_identity;
    M_dEFk = M_R.transpose()* dE * M_R;

    for(int idim = 0; idim < 3; idim++)
        for(int jdim = 0; jdim < 3; jdim++)
            M_AdEFk(idim,jdim) = M_A(idim,jdim)*M_dEFk(idim,jdim);

//    auto Jac = 2.0*M_AEFk.contract(M_dEFk) * M_Sk + M_CexpQ * M_R * M_AdEFk * M_R.transpose();
//
//    auto aux1 = Jac.contract(M_Ck) / 3.0 * M_Cinvk;
//
//    Jac = Jac + M_Sk.contract(M_Ck) / 3.0 * CinvdCCinv - aux1 - M_Sk.contract(dC) / 3.0 * M_Cinvk;
    ////    auto CinvdCCinv = M_Cinvk * dC * M_Cinvk;
//    auto Jac = 0.0 * M_identity;
//    double tmp = 0.0;
    auto Jac = 2.0*M_AEFk.contract(M_dEFk) * M_Sk + M_CexpQ * M_R * M_AdEFk * M_R.transpose();

    if(M_active)
    {
        std::cout << "Active Guccione Not Yet Coded" << std::endl;
        throw std::runtime_error("error");   }
    else
    {

    }
	/*
    *  dS4 = ( J p )' * dJ/dC * dC
    *             -   J p C^-1 dC C^-1
    *         =  0.5 * J * ( C^-1 : dC ) ( p  ) * C^-1
    *             -   J p C^-1 dC C^-1
    */
    double p = M_pressure;
    Jac += 0.5 *  p * M_Jk * M_Cinvk.contract(dC) * M_Cinvk;
    Jac -= M_Jk * p  * CinvdCCinv;
    M_deviatoric_jacobian = dU * M_total_stress + M_Fk * Jac;
}

void
Guccione::evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q)
{
    auto  dF = dU;
    auto  dC = M_Fk.transpose() * dF + dF.transpose() * M_Fk;
    auto  dE = 0.5 * dC;
    M_dEFk = M_R.transpose()* dE * M_R;

    for(int idim = 0; idim < 3; idim++)
        for(int jdim = 0; jdim < 3; jdim++)
            M_AdEFk(idim,jdim) = M_A(idim,jdim)*M_dEFk(idim,jdim);

    auto Jac = 2.0*M_AEFk.contract(M_dEFk) * M_Sk + M_CexpQ * M_R * M_AdEFk * M_R.transpose();
////    auto CinvdCCinv = M_Cinvk * dC * M_Cinvk;
//    auto Jac = 0.0 * M_identity;
//    double tmp = 0.0;
    if(M_active)
    {
        std::cout << "Active Guccione Not Yet Coded" << std::endl;
        throw std::runtime_error("error");   }
    else
    {

    }
    double p = evaluatePressure();
    double dJp = p + M_Jk * d2U(M_Jk);
    Jac += 0.5 * dJp * M_Jk * M_Cinvk.contract(dC) * M_Cinvk;
    Jac -= M_Jk * p  * M_Cinvk * dC * M_Cinvk;
    M_total_jacobian = dU * M_total_stress + M_Fk * Jac;
}


double
Guccione::evaluatePressure()
{
    double kappa = M_kappa;
    M_Fk = M_identity + M_gradU;
    M_Jk = M_Fk.det();
    //double dU = std::log(M_Jk);
    double dU = M_Jk - 1;
    return kappa * dU;
}


double
Guccione::d2U( double J)
{
    // For the time being this material is only incompressible
    // we return 1 to be compatible with the code
//    return M_kappa / M_Jk;
    return M_kappa;
}
double
Guccione::d3U( double J)
{
    return 0.0;
  //  return - M_kappa / M_Jk  / M_Jk;
}


/////// MATERIALS:
/////// CONSTITUTIVE LAWS:





} /* namespace BeatIt */



