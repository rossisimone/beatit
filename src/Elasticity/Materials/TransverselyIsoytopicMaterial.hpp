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
 * \file TransverselyIsoytopicMaterial.hpp
 *
 * \class TransverselyIsoytopicMaterial
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

#ifndef SRC_ELASTICITY_MATERIALS_TRANSVERSELYISOYTOPICMATERIAL_HPP_
#define SRC_ELASTICITY_MATERIALS_TRANSVERSELYISOYTOPICMATERIAL_HPP_

#include "Elasticity/Materials/IsotropicMaterial.hpp"

namespace BeatIt
{

enum class HexagonalCL
{
    // Transversely isotropic biological, soft tissue must be modelled using both anisotropic invariants
    // J.G. Murphy - eq (37)
    // W = mu * (I_1 - 3)  : Isotropic Part
    //   + mu5 / 2 ( 2 * I4 - I5 - 1 ) + mu4 / 2 ( I4 - 1 )^2
    // mu5 = (mu - muL)
    // mu4 = (EL + mu - 4 * muL) / 4
    // W4 = mu5 + mu4 * (I4 -1 )
    // W5 = - mu5 /2
    // W44 = mu4
    // W45 = W54 = W55 = 0
    Murphy,
	// Several law, see e.g. Holzapfel ( 2000 )
	// A new constitutive framework for arterial wall mechanics and a comparative study of material models
    // W = mu * (I_1 - 3)  : Isotropic Part
    //   + af / ( 2 bf )  ( exp( bf *  ( I4f - 1 )^2 ) - 1 )
    // W4 =  af * ( I4f - 1 ) * exp( bf *  ( I4f - 1 )^2 )
	// W5 = 0
    // W44 =  (1 + 2 * bf * ( I4f - 1 )^2 ) * af * exp( bf *  ( I4f - 1 )^2 )
	// W45 = W54 = W55 = 0
	Exponential,
	// Valve model of my invention
	// W = a1 / b1 * exp( b1 * ( I1 - 3 ) )
	//   + a2 / b2 * exp( b2 * ( I2 - 3 ) )
    //   + a4 / b4 * exp( b4 * ( I4 - 1 )^2 )
    //   + a5 / b5 * exp( b5 * ( I5 - 1 )^2 )
	// In which we have
	// W1 = a1 *  exp( b1 * ( I1 - 3 ) )
    // W11= a1 * b1 * exp( b1 * ( I1 - 3 ) )
    // W2 = a2 *  exp( b2 * ( I2 - 3 ) )
    // W22= a2 * b2 * exp( b2 * ( I2 - 3 ) )
	// W12= 0
    // W4 = a4 * (I4-1)  exp( b4 * ( I4 - 1 )^2 )
    // W5 = a5 * (I5-1)  exp( b5 * ( I5 - 1 )^2 )
    // W44= a4 * b4 * (I4-1)^2  exp( b4 * ( I4 - 1 )^2 )
    // W55= a5 * b5 * (I5-1)^2  exp( b5 * ( I5 - 1 )^2 )
	Valve
};
/*!
 *
 */
class TransverselyIsoytopicMaterial : public virtual IsotropicMaterial
{
public:
    TransverselyIsoytopicMaterial();
    ~TransverselyIsoytopicMaterial();

    void setup(GetPot& data, std::string section);
    /// This method is used for primal elasticity
    void evaluateVolumetricStress();
    void evaluateDeviatoricStress();
    void evaluateStress( ElasticSolverType solverType);
    /// This method is used only for mixed elasticy
    void evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q = 0.0) ;

    /// This method is used for primal elasticity
    void evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q = 0.0);


    double W4 (double I4, double I5, double J);
    double W44(double I4, double I5, double J);
    double W45(double I4, double I5, double J);
    double W5 (double I4, double I5, double J);
    double W54(double I4, double I5, double J);
    double W55(double I4, double I5, double J);

    HexagonalCL M_hexagonalCL;
    double M_I4;
    double M_I5;
    double M_W4;
    double M_W5;
    double M_W44;
    double M_W54;
    double M_W45;
    double M_W55;
    libMesh::VectorValue <double> M_fA;
    libMesh::TensorValue <double> M_fof;
    libMesh::TensorValue <double> M_Cfof;

    virtual void updateVariables();


};

Material* createTransverselyIsoytopicMaterial();

namespace
{
    static bool register_transversely_isotropic_material = BeatIt::Material::MaterialFactory::Register("hexagonal", &createTransverselyIsoytopicMaterial);
}


} /* namespace BeatIt */

#endif /* SRC_ELASTICITY_MATERIALS_TRANSVERSELYISOYTOPICMATERIAL_HPP_ */
