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
	M_parameters[0] = data(section+"/rho", 0.0); // rho
    M_active = data(section+"/active", false);

    double nu = data(section+"/nu", -1.0);
    if(nu <= -1 || nu >= 0.5) M_isIncompressible = true;
    else M_isIncompressible = false;

    M_C = data(section+"/C", 0.0);
    M_bff = data(section+"/bff", 0.0);
    M_bfs = data(section+"/bfs", 0.0);
    M_bt = data(section+"/bt", 0.0);

    //double nu = data(section+"/nu", 0.0);
    M_parameters[1] = E / ( 2.0 * ( 1 + nu ) );// mu
    M_parameters[2] = 1.0;// kappa
    M_kappa = M_parameters[2];
    std::cout << "\t density = " << M_parameters[0] << std::endl;
    std::cout << "\t shear modulus = " << M_parameters[1] << std::endl;
    std::cout << "\t bulk modulus = " << M_parameters[2] << std::endl;
    M_tau = 0.5 /  M_parameters[1];

    M_density = M_parameters[0];

    // stabilization
    M_tau = 0.5 / M_C;
    double tau_coeff = data(section+"/ctau", 1.0); // stabilization coefficient
    M_tau *= tau_coeff; // rho

}


void
Guccione::evaluateDeviatoricStress()
{
    /*

     */
    M_Q = Q(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);
    M_expQ = std::exp(Q);
    M_CexpQ = M_C * M_expQ;

    // Derivatives
    M_Q4f = Q4f(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);
    M_Q4s = 4s(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);
    M_Q4n = 4n(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);
    M_Q8fs = 8fs(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);
    M_Q8fn = 8fn(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);
    M_Q8sn = 8sn(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);

    //Second derivatives Q4f
    M_Q4f4f = 4f4f(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);
    //Second derivatives Q4s
    M_Q4s4s = 4s4s(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);
    //Second derivatives Q4n
    M_Q4n4n = 4n4n(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);
    //Second derivatives Q8fs
    M_Q8fs8fs = 8fs8fs(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);
    //Second derivatives Q8fn
    M_Q8fn8fn = 8fn8fn(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);
    //Second derivatives Q8sn
    M_Q8sn8sn = 8sn8sn(M_I4f, M_I4s, M_I4n, M_I8fs, M_I8fn, M_I8sn, M_Jk);

    dQ[0] = M_Q4f;
    dQ[1] = M_Q4s;
    dQ[2] = M_Q4n;
    dQ[3] = M_Q8fs;
    dQ[4] = M_Q8fn;
    dQ[5] = M_Q8sn;
    d2Q[0] = M_Q4f4f;
    d2Q[1] = M_Q4s4s;
    d2Q[2] = M_Q4n4n;
    d2Q[3] = M_Q8fs8fs;
    d2Q[4] = M_Q8fn8fn;
    d2Q[5] = M_Q8sn8sn;

    // values for loops:
    if(M_active)
    {
        std::cout << "Active Guccione Not Yet Coded" << std::endl;
        throw std::runtime_error("error");
    }
    else
    {
        // S = 0.5 * C e^Q * dQdI * dI/dE )
        M_deviatoric_stress *= 0.0;
        M_S*= 0.0;
        for(int i = 0; i < 6; i++)
        {
            M_S += 0.5 * M_CexpQ * dQ[i] * dIdE[i];
            M_deviatoric_stress += 0.5 * M_CexpQ * dQ[i] * dIdEDEV[i];
        }
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

	double p = kappa * ( M_Jk - 1);
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
            for(int jdim = 0; jdim < 3; jdim++)
            {
                //f, s, n
                M_fof(idim, jdim)  = M_f0(idim) * M_f0(jdim);
                M_sos(idim, jdim)  = M_s0(idim) * M_s0(jdim);
                M_non(idim, jdim)  = M_n0(idim) * M_n0(jdim);
                // fs
                M_fos(idim, jdim)  = M_f0(idim) * M_s0(jdim);
                M_sof(idim, jdim)  = M_s0(idim) * M_f0(jdim);
                // fn
                M_fon(idim, jdim)  = M_f0(idim) * M_n0(jdim);
                M_nof(idim, jdim)  = M_n0(idim) * M_f0(jdim);
                // sn
                M_nos(idim, jdim)  = M_n0(idim) * M_s0(jdim);
                M_son(idim, jdim)  = M_s0(idim) * M_n0(jdim);
            }
        }
        M_I4f  = M_fof.contract(M_Ck);
        M_I4s  = M_sos.contract(M_Ck);
        M_I4n  = M_non.contract(M_Ck);
        M_I8fs = M_fos.contract(M_Ck);
        M_I8fn = M_fon.contract(M_Ck);
        M_I8sn = M_son.contract(M_Ck);

        dIdE[0] = M_fof;
        dIdE[1] = M_sos;
        dIdE[2] = M_non;
        dIdE[3] = 2.0 * ( M_fos + M_sof );
        dIdE[4] = 2.0 * ( M_fon + M_nof );
        dIdE[5] = 2.0 * ( M_nos + M_son );

        // dIdEDEV = dIdE - (dIdE: C) / 3 * Cinv
        // dIdEDEV[0] = M_fof - M_I4f / 3.0 * M_Cinvk;
        // dIdEDEV[1] = M_sos - M_I4f / 3.0 * M_Cinvk;
        // dIdEDEV[2] = M_non - M_I4f / 3.0 * M_Cinvk;
        // dIdEDEV[3] = M_fos + M_sof - 2.0 * M_I8fs / 3.0 * M_Cinvk;
        // dIdEDEV[4] = M_fon + M_nof - 2.0 * M_I8fn / 3.0 * M_Cinvk;
        // dIdEDEV[5] = M_nos + M_son - 2.0 * M_I8sn / 3.0 * M_Cinvk;

        for(int l = 0; l < 6; ++l)
        {
            dIdEDEV[l] = dIdE[l] - dIdE[l].contract(M_Ck) / 3.0 * M_Cinvk;
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
Guccione::evaluateVolumetricJacobian( const libMesh::TensorValue <double>& dU, double q)
{
	M_volumetric_jacobian = q * M_Fk * M_Cinvk;
}


void
Guccione::evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q )
{
	auto  dF = dU;
	auto  dC = M_Fk.transpose() * dF + dF.transpose() * M_Fk;
    auto  dE = 0.5 * dC;
    auto CinvdCCinv = M_Cinvk * dC * M_Cinvk;
    auto Jac = 0.0 * M_identity;

    if(M_active)
    {
        std::cout << "Active Guccione Not Yet Coded" << std::endl;
        throw std::runtime_error("error");
    }
    else
    {
        // SDEV = C / 2 * exp(Q) * dQ/dI_i * (dI_i/dE)_DEV
        // dSDEV/dE : dE

        // Part 1:
        // d ( C / 2 * exp(Q) ) /dE: dE = S:dE
        // (S:dE) * 2 / C / exp(Q) * SDEV
        double SdE = M_S.contract(dE);
        Jac += 2.0 / M_CexpQ * SdE * M_deviatoric_stress;

        // Part 2:
        // Take d / dE ( dQ/dI_i ) : dE = d2Q / dI_i dI_j * (dI_j/dE:dE)
        for (int l = 0; l < 6; ++l)
        {
            Jac += 0.5 * M_CexpQ * d2Q[l] * dIdE[l].contract(dE) * dIdEDEV[l];
        }
        // Part 3
        // d/dE ((dI_i/dE)_DEV) : dE
        // 2/3 (dI_idE:C) * CinvdCCinv - 1 / 3 * (dI_idE:dC) Cinv
        for (int l = 0; l < 6; ++l)
        {
            Jac += 0.5 * M_CexpQ * dQ[l] *
                 ( 2.0 / 3.0 * dIdE[l].contract(C) * CinvdCCinv - dIdE[l].contract(dC) / 3.0 * Cinv  );
        }
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
Guccione::evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q)
{
    auto  dF = dU;
    auto  dC = M_Fk.transpose() * dF + dF.transpose() * M_Fk;
    auto  dE = 0.5 * dC;
    auto CinvdCCinv = M_Cinvk * dC * M_Cinvk;
    auto Jac = 0.0 * M_identity;
    if(M_active)
    {
        std::cout << "Active Guccione Not Yet Coded" << std::endl;
        throw std::runtime_error("error");   }
    else
    {
        // SDEV = C / 2 * exp(Q) * dQ/dI_i * (dI_i/dE)_DEV
        // dSDEV/dE : dE

        // Part 1:
        // d ( C / 2 * exp(Q) ) /dE: dE = S:dE
        // (S:dE) * 2 / C / exp(Q) * SDEV
        double SdE = M_S.contract(dE);
        Jac += 2.0 / M_CexpQ * SdE * M_deviatoric_stress;

        // Part 2:
        // Take d / dE ( dQ/dI_i ) : dE = d2Q / dI_i dI_j * (dI_j/dE:dE)
        for (int l = 0; l < 6; ++l)
        {
            Jac += 0.5 * M_CexpQ * d2Q[l] * dIdE[l].contract(dE) * dIdEDEV[l];
        }
        // Part 3
        // d/dE ((dI_i/dE)_DEV) : dE
        // 2/3 (dI_idE:C) * CinvdCCinv - 1 / 3 * (dI_idE:dC) Cinv
        for (int l = 0; l < 6; ++l)
        {
            Jac += 0.5 * M_CexpQ * dQ[l] *
                 ( 2.0 / 3.0 * dIdE[l].contract(C) * CinvdCCinv - dIdE[l].contract(dC) / 3.0 * Cinv  );
        }
    }
    double p = evaluatePressure();
    double dJp = p + M_Jk * d2U(M_Jk);
    Jac += 0.5 * dJp * M_Jk * M_Cinvk.contract(dC) * M_Cinvk;
    Jac -= M_Jk * p  * M_Cinvk * dC * M_Cinvk;
    M_total_jacobian = dU * M_total_stress + M_Fk * Jac;
    //std::cout << "--- ISOTROPIC MATERIAL: evaluateJacobian not coded for isotropic material!!!" << std::endl;
   // throw std::runtime_error("Primal Formulation not coded for isotropic material");
}


double
Guccione::evaluatePressure()
{
    double kappa = M_kappa;
    M_Fk = M_identity + M_gradU;
    M_Jk = M_Fk.det();
    double dU = std::log(M_Jk);
    return kappa * dU;
}

double
Guccione::evaluatePressureResidual()
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
	    //std::cout << "--- ISOTROPIC MATERIAL: evaluatePressureResidual compressible case not coded for isotropic material!!!" << std::endl;
	    //throw std::runtime_error("Primal Formulation not coded for isotropic material");
	}
	return - RHS;

}

double
Guccione::dpdF(const libMesh::TensorValue <double>&  dF)
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
        //std::cout << "--- ISOTROPIC MATERIAL: dpdF compressible case not coded for isotropic material!!!" << std::endl;
        //throw std::runtime_error("Primal Formulation not coded for isotropic material");
		//return 0.0;
	}
}

double
Guccione::d2U( double J)
{
    // For the time being this material is only incompressible
    // we return 1 to be compatible with the code
    return M_kappa / M_Jk;
}
double
Guccione::d3U( double J)
{
    return - M_kappa / M_Jk  / M_Jk;
}


/////// MATERIALS:
/////// CONSTITUTIVE LAWS:





// Helper parameters
double
Guccione::Q (double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.25 * ( bff * (I4f-1) * (I4f-1)
                  + bt  * ( (I4s-1) * (I4s-1) + (I4n-1) * (I4n-1) + 2 * I8sn * I8sn  )
              + 2 * bfs * ( I8fs * I8fs + I8fn * I8fn )   );
}
// Derivatives
double Guccione::Q4f(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bff * (I4f-1);
}
double Guccione::Q4s(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bt * (I4s-1);
}

double Guccione::Q4n(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bt * (I4n-1);
}

double Guccione::Q8fs(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bfs * I8fs;
}

double Guccione::Q8fn(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bfs * I8fn;
}

double Guccione::Q8sn(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bt * I8sn;
}


//Second derivatives Q4f
double Guccione::Q4f4f(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bff;
}

//Second derivatives Q4s
double Guccione::Q4s4s(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bt;
}

//Second derivatives Q4n
double Guccione::Q4n4n(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bt;
}

//Second derivatives Q8fs
double Guccione::Q8fs8fs(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bfs;
}

//Second derivatives Q8fn
double Guccione::Q8fn8fn(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bfs;
}

//Second derivatives Q8sn
double Guccione::Q8sn8sn(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J)
{
    return 0.5 * M_bt;
}



} /* namespace BeatIt */



