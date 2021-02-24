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
 * Material.cpp
 *
 *  Created on: Oct 21, 2016
 *      Author: srossi
 */

#include "Elasticity/Materials/Material.hpp"

namespace BeatIt {

namespace MaterialUtilities
{

void cross_product(const libMesh::TensorValue <double>& A, const libMesh::TensorValue <double>& B, libMesh::TensorValue <double>& C )
{
    C(0,0) = A(1,1) * B(2,2) - A(1,2) * B(2,1) - A(2,1) * B(1,2) + A(2,2) * B(1,1);
    C(0,1) = A(1,2) * B(2,0) - A(2,2) * B(1,0) - A(1,0) * B(2,2) + A(2,0) * B(1,2);
    C(0,2) = A(1,0) * B(2,1) - A(2,0) * B(1,1) - A(1,1) * B(2,0) + A(2,2) * B(1,0);

    C(1,0) = A(2,1) * B(0,2) - A(0,1) * B(2,2) - A(2,2) * B(0,1) + A(0,2) * B(2,1);
    C(1,1) = A(2,2) * B(0,0) - A(2,0) * B(0,2) - A(0,2) * B(2,0) + A(0,0) * B(2,2);
    C(1,2) = A(2,0) * B(0,1) - A(0,0) * B(2,1) - A(2,1) * B(0,0) + A(0,1) * B(2,0);

    C(2,0) = A(0,1) * B(1,2) - A(1,1) * B(0,2) - A(0,2) * B(1,1) + A(1,2) * B(0,1);
    C(2,1) = A(0,2) * B(1,0) - A(1,2) * B(0,0) - A(0,0) * B(1,2) + A(1,0) * B(0,2);
    C(2,2) = A(0,0) * B(1,1) - A(0,1) * B(1,0) - A(1,0) * B(0,1) + A(1,1) * B(0,0);

}

void cof(const libMesh::TensorValue <double>& Fk, libMesh::TensorValue <double>& H )
{
    cross_product(Fk, Fk, H);
    H *= 0.5;
}

void dH(const libMesh::TensorValue <double>& dU, const libMesh::TensorValue <double>& Fk, libMesh::TensorValue <double>& dH )
{
    cross_product(Fk, dU, dH);
}

} // MaterialUtilities


const libMesh::TensorValue <double>  Material::M_identity
	= libMesh::TensorValue <double>( 1.0, 0.0, 0.0,
																  0.0, 1.0, 0.0,
																  0.0, 0.0, 1.0 );

Material::Material()
	: M_PK1()
	, M_total_stress()
	, M_deviatoric_stress()
	, M_volumetric_stress()
	, M_total_jacobian()
	, M_deviatoric_jacobian()
	, M_volumetric_jacobian()
	, M_gradU()
	, M_strain()
	, M_strain_rate()
	, M_Fk()
    , M_Hk()
	, M_Ck()
	, M_Cinvk()
	, M_Jk(1.0)
	, M_pressure()
	, M_parameters()
	, M_f0()
	, M_isIncompressible(false)
	, M_ndim(3)
	, M_density(1.0)
	, M_tau(0.0)
    , M_FA()
{

}

Material::~Material() {
	// TODO Auto-generated destructor stub
}

void
Material::setup(GetPot& data, std::string section, int ndim)
{
	M_ndim = ndim;
	setup(data, section);
}



void
Material::dH(const libMesh::TensorValue <double>&  dU, libMesh::TensorValue <double> dcof)   const
{
    MaterialUtilities::cross_product(M_Fk, dU, dcof);
}

void
Material::evaluateStress( ElasticSolverType solverType)
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
            M_total_stress += M_Jk * M_pressure * M_Cinvk;
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

//  double mu = M_parameters[1];
//  M_PK1 = mu * (M_Fk - M_identity) +  M_Jk * M_pressure * M_Cinvk;
}

void
Material::evaluateVolumetricJacobian( const libMesh::TensorValue <double>& dU, double q)
{
    // This term represets
    // Pdev = J p F^-T
    // J dp F^-T * phi
    M_volumetric_jacobian = M_Jk * q * M_Fk * M_Cinvk;
}

double
Material::evaluatePressureResidual()
{
    /*
     *  Evaluate the RHS of the pressure equation with a minus, that is
     *  p = RHS
     *  p - RHS = 0
     *  return -RHS
     */
    double RHS;
    if(M_isIncompressible)
    {
        RHS = ( M_Jk - 1.0);
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
Material::dpdF(const libMesh::TensorValue <double>&  dF)
{
    /*
     *   p - k U'(J) = 0
     *
     *   dp / dF = k U''(J) dJ/DF = k U''(J) J F^-T
     *
     *  dPvol/dp =
     */
    /*
     *   p - k U'(J) = 0
     *
     *   dp / dF = k U''(J) dJ/DF = k U''(J) J F^-T
     *
     *  dPvol/dp =
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


} /* namespace BeatIt */
