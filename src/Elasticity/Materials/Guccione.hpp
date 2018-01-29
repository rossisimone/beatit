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
	void evaluateVolumetricJacobian( const libMesh::TensorValue <double>& dU, double q = 0.0);

	/// This method is used for primal elasticity
	void evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q = 0.0);

	double evaluatePressure();
    double evaluatePressureResidual();
    double dpdF(const libMesh::TensorValue <double>&  dF);

    double d2U( double J = 1.0);
    double d3U( double J = 1.0);


    // Helper parameters
    double Q (double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);

    // Derivatives
    double Q4f(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);
    double Q4s(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);
    double Q4n(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);
    double Q8fs(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);
    double Q8fn(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);
    double Q8sn(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);

    //Second derivatives Q4f
    double Q4f4f(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);

    //Second derivatives Q4s
    double Q4s4s(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);

    //Second derivatives Q4n
    double Q4n4n(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);

    //Second derivatives Q8fs
    double Q8fs8fs(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);

    //Second derivatives Q8fn
    double Q8fn8fn(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);

    //Second derivatives Q8sn
    double Q8sn8sn(double I4f, double I4s, double I4s, double I8fs, double I8fn, double I8sn, double J);


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

    // Derivatives
    double M_Q4f;
    double M_Q4s;
    double M_Q4n;
    double M_Q8fs;
    double M_Q8fn;
    double M_Q8sn;


    //Second derivatives Q4f
    double M_Q4f4f;

    //Second derivatives Q4s
    double M_Q4s4s;

    //Second derivatives Q4n
    double M_Q4n4n;

    //Second derivatives Q8fs
    double M_Q8fs8fs;

    //Second derivatives Q8fn
    double M_Q8fn8fn;

    //Second derivatives Q8sn
    double M_Q8sn8sn;

    double M_kappa;
    bool M_active;

    libMesh::TensorValue <double> M_fof;
    libMesh::TensorValue <double> M_sos;
    libMesh::TensorValue <double> M_non;
    libMesh::TensorValue <double> M_fos;
    libMesh::TensorValue <double> M_sof;
    libMesh::TensorValue <double> M_fon;
    libMesh::TensorValue <double> M_nof;
    libMesh::TensorValue <double> M_nos;
    libMesh::TensorValue <double> M_son;

    libMesh::TensorValue <double> M_S;
    // ordered as f, s, n, fs, fn, sn
    double dQ[6];
    // ordered as f, s, n, fs, fn, sn
    double d2Q[6];
    // ordered as fof, sos, non, sym(fos), sym(fon), sym(son)
    libMesh::TensorValue <double> dIdE[6]; // dIdE = 2 * dIdC
    libMesh::TensorValue <double> dIdEDEV[6]; // dIdE = dIdE - 1 / (3J) (dIdE:C) Cinv

    virtual void updateVariables();

};



Material* createGuccione();

namespace
{
	static bool register_isotropic_material = BeatIt::Material::MaterialFactory::Register("guccione", &createGuccione);
}


} /* namespace BeatIt */

#endif /* SRC_ELASTICITY_MATERIALS_Guccione_HPP_ */
