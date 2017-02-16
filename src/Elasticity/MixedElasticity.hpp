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
 * MixedElasticity.hpp
 *
 *  Created on: Oct 24, 2016
 *      Author: srossi
 */

#ifndef SRC_ELASTICITY_MIXEDELASTICITY_HPP_
#define SRC_ELASTICITY_MIXEDELASTICITY_HPP_

#include "Elasticity/Elasticity.hpp"

namespace BeatIt {

class MixedElasticity : public virtual Elasticity
{
public:
	MixedElasticity( libMesh::EquationSystems& es, std::string system_name );
	virtual ~MixedElasticity();

    void assemble_residual(double dt = 0.0, libMesh::NumericVector<libMesh::Number>* activation_ptr = nullptr);
    //void assemble_residual(libMesh::NumericVector<libMesh::Number>* activation_ptr = nullptr);
    void assemble_jacobian(){}
};

} /* namespace BeatIt */

#endif /* SRC_ELASTICITY_MIXEDELASTICITY_HPP_ */
