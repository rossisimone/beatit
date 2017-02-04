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
	// TODO Auto-generated constructor stub

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


} /* namespace BeatIt */
