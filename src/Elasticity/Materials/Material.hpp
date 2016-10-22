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
#include "Util/Factory.hpp"
#include <vector>


/// Forward Declarations
class GetPot;

namespace BeatIt {

class Material {
public:

	typedef Factory<Material, std::string>     MaterialFactory;
	Material();
	virtual ~Material();

	virtual void setup(GetPot& data, std::string section) = 0;

	virtual void evaluateVolumetricStress() = 0;
	virtual void evaluateDeviatoricStress() = 0;
	virtual void evaluateStress() = 0;

	virtual void evaluateVolumetricJacobian( const libMesh::TensorValue <double>& dU, double q = 0.0) = 0;
	virtual void evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q = 0.0) = 0;
	virtual void evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q = 0.0) = 0;


   libMesh::TensorValue <double> M_total_stress;
   libMesh::TensorValue <double> M_deviatoric_stress;
   libMesh::TensorValue <double> M_volumetric_stress;
   libMesh::TensorValue <double> M_total_jacobian;
   libMesh::TensorValue <double> M_deviatoric_jacobian;
   libMesh::TensorValue <double> M_volumetric_jacobian;
   libMesh::TensorValue <double> M_gradU;
   libMesh::TensorValue <double> M_strain;
   libMesh::TensorValue <double> M_strain_rate;
   libMesh::TensorValue <double> M_Fk;
   libMesh::TensorValue <double> M_Ck;
   libMesh::TensorValue <double> M_Cinvk;
   double M_Jk;
   double M_pressure;
   std::vector<double>  M_parameters;
   std::vector<double>  M_fibers;

   const static libMesh::TensorValue <double> M_identity;

};

} /* namespace BeatIt */

#endif /* SRC_ELASTICITY_MATERIALS_MATERIAL_HPP_ */
