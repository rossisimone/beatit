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

/**
 * \file HolzapfelOgden.cpp
 *
 * \class HolzapfelOgden
 *
 * \brief This class provides a simple factory implementation
 *
 * For details on how to use it check the test_factory in the testsuite folder
 *
 *
 * \author srossi
 *
 * \version 0.0
 *
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Dec 29, 2016
 *
 */

#include "Elasticity/Materials/HolzapfelOgden.hpp"
#include "libmesh/getpot.h"

namespace BeatIt
{


Material* createHolzapfelOgden()
{
    return new HolzapfelOgden;
}


HolzapfelOgden::HolzapfelOgden()
 : M_I4f(1.0)
 , M_I4s(1.0)
 , M_I8fs(0.0)
 , M_W4f(0.0)
 , M_W44f(0.0)
 , M_W4s(0.0)
 , M_W44s(0.0)
 , M_W8fs(0.0)
 , M_W88fs(0.0)
{
    /*
     *  The parameters will be stored using the following order
     *  0: density = rho
     *  1: shear modulus: a
     *  2: bulk modulus: b
     *  3: shear modulus: a4f
     *  4: bulk modulus: b4f
     *  5: shear modulus: a4s
     *  6: bulk modulus: b4s
     *  7: shear modulus: a8fs
     *  8: bulk modulus: b8fs
     */
    M_parameters.resize(9);
}

HolzapfelOgden::~HolzapfelOgden()
{
    // TODO Auto-generated destructor stub
}


void
HolzapfelOgden::setup(GetPot& data, std::string section)
{
    IsotropicMaterial::setup(data, section);
	// From the isotropic material
	// M_parameters[0] = density;
	// M_parameters[1] = a;
	// M_parameters[2] = b;

	double a4f = data(section+"/af", 0.0);
	double b4f = data(section+"/bf", 0.0);
	M_parameters[3] = a4f;
	M_parameters[4] = b4f;

	double a4s = data(section+"/as", 0.0);
	double b4s = data(section+"/bs", 0.0);
	M_parameters[5] = a4s;
	M_parameters[6] = b4s;

	double a8fs = data(section+"/afs", 0.0);
	double b8fs = data(section+"/bfs", 0.0);
	M_parameters[7] = a8fs;
	M_parameters[8] = b8fs;
}

/// This method is used for primal elasticity
void
HolzapfelOgden::evaluateVolumetricStress()
{
 // NOT IMPLEMENTED
 // ONLY INCOMPRESSIBLE
    std::cout << "--- HolzapfelOgden : Primal Formulation not coded for HolzapfelOgden !!!" << std::endl;
    throw std::runtime_error("Primal Formulation not coded for  HolzapfelOgden");

}
void
HolzapfelOgden::evaluateDeviatoricStress()
{
    IsotropicMaterial::evaluateDeviatoricStress();

    double a4f = 	M_parameters[3];
	double b4f = M_parameters[4];
	double a4s =	M_parameters[5];
	double b4s = M_parameters[6];
	double a8fs =M_parameters[7];
	double b8fs = M_parameters[8];
    M_W4f  = W4 (a4f, b4f, M_I4f, M_Jk);
    M_W44f = W44(a4f, b4f, M_I4f, M_Jk);
    M_W4s  = W4 (a4s, b4s, M_I4s, M_Jk);
    M_W44s = W44(a4s, b4s, M_I4s, M_Jk);
    M_W8fs  = W8 (a8fs, b8fs, M_I8fs, M_Jk);
    M_W88fs = W88(a8fs, b8fs, M_I8fs, M_Jk);
    //std::cout << "Stress: " << M_W4f << std::endl;
    // Fiber part
    //M_deviatoric_stress += 2.0 * M_W4f * M_fof;
    M_deviatoric_stress += 2.0 * M_W4f * ( M_fof - M_I4f / 3.0 * M_Cinvk);
    //M_deviatoric_stress += 2.0 * M_W4f * ( M_fof - M_I4f / 3.0 * M_Cinvk);
	// Sheet part
    M_deviatoric_stress += 2.0 * M_W4s * ( M_sos - M_I4s / 3.0 * M_Cinvk);
    //M_deviatoric_stress += 2.0 * M_W4s * ( M_sos - M_I4s / 3.0 * M_Cinvk);
	// Sheet part
    //M_deviatoric_stress += M_W8fs * ( M_fos + M_sof );
    M_deviatoric_stress += M_W8fs * ( M_fos + M_sof - 2.0 * M_I8fs / 3.0 * M_Cinvk);


}


void
HolzapfelOgden::updateVariables()
{
    IsotropicMaterial::updateVariables();
    // f = F f0
    // but in order to compute I4 and I5 without defining different quantities I use f = C f0
    if(M_active)
    {
        std::cout << "Active Holzapfel Ogden Not Yet Coded" << std::endl;
        throw std::runtime_error("error");
//        M_fA  = M_FAinv * M_f0;
//        for(int idim = 0; idim < 3; idim++)
//        {
//            for(int jdim = 0; jdim < 3; jdim++)
//            {
//                M_fof(idim, jdim)  = M_fA(idim) * M_fA(jdim);
//            }
//
//        }
//        M_I4 = M_fof.contract(M_Ck);
//        M_I5 = M_fof.contract(M_Ck * M_CAinv * M_Ck);
//        M_Cfof = M_CAinv * M_Ck * M_fof;
    }
    else
    {
        for(int idim = 0; idim < 3; idim++)
        {
            for(int jdim = 0; jdim < 3; jdim++)
            {
                M_fof(idim, jdim)  = M_f0(idim) * M_f0(jdim);
                M_fos(idim, jdim)  = M_f0(idim) * M_s0(jdim);
                M_sof(idim, jdim)  = M_s0(idim) * M_f0(jdim);
                M_sos(idim, jdim)  = M_s0(idim) * M_s0(jdim);
            }

        }
        M_I4f   = M_fof.contract(M_Ck);
        M_I4s  = M_sos.contract(M_Ck);
		M_I8fs = M_fos.contract(M_Ck);
    }
}




/// This method is used only for mixed elasticy
void
HolzapfelOgden::evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q )
{
    IsotropicMaterial::evaluateDeviatoricJacobian(dU, q);
    auto  dF = dU;
    auto  dC = M_Fk.transpose() * dF + dF.transpose() * M_Fk;
    auto CinvdCCinv = M_Cinvk * dC * M_Cinvk;
    auto Jac = 0.0 * M_identity;

    // I4f part
     auto dI4f = dC.contract(M_fof);
     auto dW4f = M_W44f * dI4f;
     Jac += 2.0 * dW4f * ( M_fof - M_I4f / 3.0 * M_Cinvk );
     Jac += 2.0 * M_W4f * ( CinvdCCinv * M_I4f / 3.0 - dI4f / 3.0 * M_Cinvk );

     // I4s part
      auto dI4s = dC.contract(M_sos);
      auto dW4s = M_W44s * dI4s;
      Jac += 2.0 * dW4s * ( M_sos - M_I4s / 3.0 * M_Cinvk );
      Jac += 2.0 * M_W4s * ( CinvdCCinv * M_I4s / 3.0 - dI4s / 3.0 * M_Cinvk );

      // I8fs part
       auto dI8fs = dC.contract(M_fos);
       auto dW8fs = M_W88fs * dI8fs;
       Jac += 2.0 * dW8fs * ( M_fos + M_sof - 2.0 * M_I8fs / 3.0 * M_Cinvk);
       Jac += 2.0 * M_W8fs * ( 2.0 * CinvdCCinv * M_I8fs / 3.0 - 2.0 * dI8fs / 3.0 * M_Cinvk );

     M_deviatoric_jacobian += M_Fk * Jac;
}


/// This method is used for primal elasticity
void
HolzapfelOgden::evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q)
{
    // NOT IMPLEMENTED
    // ONLY INCOMPRESSIBLE
       std::cout << "--- HolzapfelOgden : evaluateJacobian not coded for HolzapfelOgden !!!" << std::endl;
       throw std::runtime_error("evaluateJacobian not coded for  HolzapfelOgden");
}


double
HolzapfelOgden::W4 (double a, double b, double I4, double J)
{
    double W4 = 0.0;
    if(I4 < 0.999 ) return 0.0;
	W4 = a * (I4 - 1) * std::exp( b * (I4-1) * (I4-1) );
    return W4;
}

double
HolzapfelOgden::W44(double a, double b, double I4, double J)
{
    double W44 = 0.0;
    if(I4 < 0.999 ) return 0.0;
	double aux = b * (I4-1) * (I4-1);
	W44 = a * std::exp( aux ) * ( 1.0 + 2.0 * aux );
    return W44;
}


double
HolzapfelOgden::W8 (double a, double b, double I8, double J)
{
    double W8 = 0.0;
	W8 = a * I8 * std::exp( b * I8 * I8 );
    return W8;
}

double
HolzapfelOgden::W88(double a, double b, double I8, double J)
{
    double W88 = 0.0;
	double aux = ( b * I8 * I8 );
	W88 = a * std::exp( aux ) * ( 1.0 + 2.0 * aux );
    return W88;
}

} /* namespace BeatIt */
