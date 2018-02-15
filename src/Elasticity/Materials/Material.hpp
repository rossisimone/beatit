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
 * Material.hpp
 *
 *  Created on: Oct 21, 2016
 *      Author: srossi
 */


#ifndef SRC_ELASTICITY_MATERIALS_MATERIAL_HPP_
#define SRC_ELASTICITY_MATERIALS_MATERIAL_HPP_

#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "Util/Factory.hpp"
#include <vector>
#include "Elasticity/ElasticSolverType.hpp"


/// Forward Declarations
class GetPot;


namespace BeatIt {

namespace MaterialUtilities
{
void cross_product(const libMesh::TensorValue <double>& A, const libMesh::TensorValue <double>& B, libMesh::TensorValue <double>& C );
void cof(const libMesh::TensorValue <double>& Fk, libMesh::TensorValue <double>& H );
void dH(const libMesh::TensorValue <double>& dU, const libMesh::TensorValue <double>& Fk, libMesh::TensorValue <double>& dH );
} // MaterialUtilities


class Material {
public:

	typedef Factory<Material, std::string>     MaterialFactory;
	Material();
	virtual ~Material();

	void setup(GetPot& data, std::string section, int ndim);
	virtual void setup(GetPot& data, std::string section) = 0;



	virtual void evaluateVolumetricStress() = 0;
	virtual void evaluateDeviatoricStress() = 0;
	virtual void evaluateStress( ElasticSolverType solverType);
	// Evaluate First Piola-Kirchhoff Stress Tensor from PK2
	// Paperinik 1


	/// This method is used only for mixed elasticy
    virtual void evaluateVolumetricJacobian( const libMesh::TensorValue <double>& dU, double q = 0.0);
	/// This method is used only for mixed elasticy
	virtual void evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q = 0.0) = 0;
	/// This method is used for primal elasticity
	virtual void evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q = 0.0) = 0;

	virtual double evaluatePressure() = 0;
	virtual double evaluatePressureResidual();
	virtual double dpdF(const libMesh::TensorValue <double>&  dF);

	//virtual double dU ( double J) = 0;
    virtual double d2U( double J = 1.0) = 0;
    virtual double d3U( double J = 1.0) = 0;
    libMesh::TensorValue <double> H()  const { return M_Hk; }
    virtual void dH(const libMesh::TensorValue <double>&  dU, libMesh::TensorValue <double> dcof)  const;

	bool isIncompressible ()
	{
		return M_isIncompressible;
	}

   libMesh::TensorValue <double> M_PK1;                         // First Piola-Kirchhoff Stress
   libMesh::TensorValue <double> M_total_stress;	      //  Second Piola-Kirchhoff Stress
   libMesh::TensorValue <double> M_deviatoric_stress;//  Second Deviatoric Piola-Kirchhoff Stress
   libMesh::TensorValue <double> M_volumetric_stress;// Second Volumetric Piola-Kirchhoff Stress
   libMesh::TensorValue <double> M_total_jacobian;       // Total Jacobian
   libMesh::TensorValue <double> M_deviatoric_jacobian;
   libMesh::TensorValue <double> M_volumetric_jacobian;
   libMesh::TensorValue <double> M_gradU;
   libMesh::TensorValue <double> M_strain;
   libMesh::TensorValue <double> M_strain_rate;
   libMesh::TensorValue <double> M_Fk;
   libMesh::TensorValue <double> M_Hk;   // Cofactor of Fk
   libMesh::TensorValue <double> M_Ck;
   libMesh::TensorValue <double> M_Cinvk;
   double M_Jk;
   double M_pressure;
   double M_pressureResidual;
   std::vector<double>  M_parameters;
   libMesh::VectorValue <double> M_f0;
   libMesh::VectorValue <double> M_s0;
   libMesh::VectorValue <double> M_n0;
   bool M_isIncompressible;
   int M_ndim;
   libMesh::TensorValue <double> M_FA;
   libMesh::TensorValue <double> M_FAinv;
   libMesh::TensorValue <double> M_CAinv;

   double M_density;
   double M_tau;

   const static libMesh::TensorValue <double> M_identity;

   virtual void updateVariables() {}

};

} /* namespace BeatIt */

#endif /* SRC_ELASTICITY_MATERIALS_MATERIAL_HPP_ */
