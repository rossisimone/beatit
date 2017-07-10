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
 * \file HolzapfelOgden.hpp
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

#ifndef SRC_ELASTICITY_MATERIALS_HolzapfelOgden_HPP_
#define SRC_ELASTICITY_MATERIALS_HolzapfelOgden_HPP_

#include "Elasticity/Materials/IsotropicMaterial.hpp"

namespace BeatIt
{

/*!
 *
 */
class HolzapfelOgden : public virtual IsotropicMaterial
{
public:
    HolzapfelOgden();
    ~HolzapfelOgden();

    void setup(GetPot& data, std::string section);
    /// This method is used for primal elasticity
    void evaluateVolumetricStress();
    void evaluateDeviatoricStress();
    void evaluateStress( ElasticSolverType solverType);
    /// This method is used only for mixed elasticy
    void evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q = 0.0) ;

    /// This method is used for primal elasticity
    void evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q = 0.0);


    double W4   (double a, double b, double I4  , double J);
    double W44(double a, double b, double I4  , double J);
    double W8  (double a, double b, double I8,  double J);
    double W88(double a, double b, double I8, double J);

    double M_I4f;
    double M_I4s;
    double M_I8fs;
    double M_W4f;
    double M_W44f;
    double M_W4s;
    double M_W44s;
	double M_W8fs;
	double M_W88fs;

    libMesh::VectorValue <double> M_fA;
    libMesh::TensorValue <double> M_fof;
    libMesh::TensorValue <double> M_sos;
    libMesh::TensorValue <double> M_fos;
    libMesh::TensorValue <double> M_sof;

    virtual void updateVariables();


};

Material* createHolzapfelOgden();

namespace
{
    static bool register_ho_material = BeatIt::Material::MaterialFactory::Register("HO", &createHolzapfelOgden);
}


} /* namespace BeatIt */

#endif /* SRC_ELASTICITY_MATERIALS_TRANSVERSELYISOYTOPICMATERIAL_HPP_ */
