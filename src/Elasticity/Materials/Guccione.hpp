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
 * Guccione.hpp
 *
 *  Created on: Oct 21, 2016
 *      Author: srossi
 */

#ifndef SRC_ELASTICITY_MATERIALS_Guccione_HPP_
#define SRC_ELASTICITY_MATERIALS_Guccione_HPP_

#include "Elasticity/Materials/Material.hpp"

namespace BeatIt {

// This class is specifically for
//
// W = C/2 (e^Q-1), Q = bf  * Eff^2
//                    + bt  * (Ess^2 + Enn^2 + Esn^2 + Ens^2 )
//                    + bfs * (Efs^2 + Esf^2 + Efn^2 + Enf^2)
//
// which we rewrite as
//
// W = C/2 (e^Q-1), Q =0.25 * ( bf  * (I4f-1)^2
//                            + bt  * ( (I4s-1)^2 + (I4n-1)^2 + 2 I8sn^2 )
//                            + bfs * ( 2 I8fs^2  + 2 I8fn^2) )

class Guccione : public virtual Material
{
public:
	Guccione();
	virtual ~Guccione();


    void setup(GetPot& data, std::string section);

	/// This method is used for primal elasticity
	void evaluateVolumetricStress();
	void evaluateDeviatoricStress();
	void evaluateStress( ElasticSolverType solverType);
	/// This method is used only for mixed elasticy
	void evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q = 0.0) ;

	/// This method is used for primal elasticity
	void evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q = 0.0);

	double evaluatePressure();

    double d2U( double J = 1.0);
    double d3U( double J = 1.0);





    // Invariants
    double M_I4f;
    double M_I4s;
    double M_I4n;
    double M_I8fs;
    double M_I8fn;
    double M_I8sn;
    double M_I3;
    double M_J2;

    // Parameters
    double M_C;
    double M_bff;
    double M_bfs;
    double M_bt;

    // Helper parameters
    double M_Q;
    double M_expQ;
    double M_CexpQ;

    double M_kappa;
    bool M_active;

    libMesh::TensorValue <double> M_Sk;
    libMesh::TensorValue <double> M_R; // Rotation Matrix
    libMesh::TensorValue <double> M_Ek; // Rotation Matrix
    libMesh::TensorValue <double> M_EFk; // Rotation Matrix
    libMesh::TensorValue <double> M_dEFk; // Rotation Matrix
    libMesh::TensorValue <double> M_A; // Rotation Matrix
    libMesh::TensorValue <double> M_AEFk; // Rotation Matrix
    libMesh::TensorValue <double> M_AdEFk; // Rotation Matrix

    virtual void updateVariables();

};



Material* createGuccione();

namespace
{
	static bool register_guccione_material = BeatIt::Material::MaterialFactory::Register("guccione", &createGuccione);
}


} /* namespace BeatIt */

#endif /* SRC_ELASTICITY_MATERIALS_Guccione_HPP_ */
