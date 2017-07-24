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
 * \file TransverselyIsoytopicMaterial.cpp
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

#include "Elasticity/Materials/TransverselyIsoytopicMaterial.hpp"
#include "libmesh/getpot.h"

namespace BeatIt
{



Material* createTransverselyIsoytopicMaterial()
{
    return new TransverselyIsoytopicMaterial;
}


TransverselyIsoytopicMaterial::TransverselyIsoytopicMaterial()
{
    /*
     *  The parameters will be stored using the following order
     *  0: density = rho
     *  1: shear modulus: mu
     *  2: bulk modulus: kappa
     */
    M_parameters.resize(5);
}

TransverselyIsoytopicMaterial::~TransverselyIsoytopicMaterial()
{
    // TODO Auto-generated destructor stub
}


void
TransverselyIsoytopicMaterial::setup(GetPot& data, std::string section)
{
    IsotropicMaterial::setup(data, section);

    std::string type = data(section+"/hexagonal", "UNKNOWN");
    std::map<std::string, HexagonalCL> hexagonal_cl_map;
    hexagonal_cl_map["murphy"] = HexagonalCL::Murphy;
    hexagonal_cl_map["exponential"] = HexagonalCL::Exponential;
    auto it_cl = hexagonal_cl_map.find(type);
    if( it_cl != hexagonal_cl_map.end() ) M_hexagonalCL = it_cl->second;
    else
    {
        std::cout << "--- TransverselyIsoytopicMaterial: hexagonal type " << type << " is not known" << std::endl;
        std::cout << "    available hexagonal types: " << std::endl;
        for(auto& j : hexagonal_cl_map) std::cout << "     " << j.first << std::endl;
        throw std::runtime_error("TransverselyIsoytopicMaterial setup error");
    }

    switch(M_hexagonalCL)
    {
        case HexagonalCL::Murphy:
        {
            double mu4 = data(section+"/mu4", 0.0);
            double mu5 = data(section+"/mu5", 0.0);
            M_parameters[3] = mu4;
            M_parameters[4] = mu5;
            break;
        }
        default:
        {
        case HexagonalCL::Exponential:
        {
            double af = data(section+"/af", 0.0);
            double bf = data(section+"/bf", 0.0);
            M_parameters[3] = af;
            M_parameters[4] = bf;
            break;
        }
            break;
        }
    }

}
/// This method is used for primal elasticity
void
TransverselyIsoytopicMaterial::evaluateVolumetricStress()
{
 // NOT IMPLEMENTED
 // ONLY INCOMPRESSIBLE
    std::cout << "--- TransverselyIsoytopicMaterial : Primal Formulation not coded for TransverselyIsoytopicMaterial !!!" << std::endl;
    throw std::runtime_error("Primal Formulation not coded for  TransverselyIsoytopicMaterial");

}
void
TransverselyIsoytopicMaterial::evaluateDeviatoricStress()
{
    IsotropicMaterial::evaluateDeviatoricStress();
    M_W4  = W4 (M_I4, M_I5, M_Jk);
    M_W44 = W44(M_I4, M_I5, M_Jk);
    M_W45 = W45(M_I4, M_I5, M_Jk);
    M_W5  = W5 (M_I4, M_I5, M_Jk);
    M_W54 = W54(M_I4, M_I5, M_Jk);
    M_W55 = W55(M_I4, M_I5, M_Jk);

    M_deviatoric_stress += 2.0 * M_W4 * ( M_fof - M_I4 / 3.0 * M_Cinvk);
    M_deviatoric_stress += 2.0 * M_W5 * ( M_Cfof + M_Cfof.transpose() - 2.0 * M_I5 / 3.0 * M_Cinvk);


}


void
TransverselyIsoytopicMaterial::updateVariables()
{
    IsotropicMaterial::updateVariables();
    // f = F f0
    // but in order to compute I4 and I5 without defining different quantities I use f = C f0
    if(M_active)
    {
//        std::cout << "WRONG PLACE" << std::endl;
//        throw std::runtime_error("error");
        M_fA  = M_FAinv * M_f0;
        for(int idim = 0; idim < 3; idim++)
        {
            for(int jdim = 0; jdim < 3; jdim++)
            {
                M_fof(idim, jdim)  = M_fA(idim) * M_fA(jdim);
            }

        }
        M_I4 = M_fof.contract(M_Ck);
        M_I5 = M_fof.contract(M_Ck * M_CAinv * M_Ck);
        M_Cfof = M_CAinv * M_Ck * M_fof;
    }
    else
    {
        for(int idim = 0; idim < 3; idim++)
        {
            for(int jdim = 0; jdim < 3; jdim++)
            {
                M_fof(idim, jdim)  = M_f0(idim) * M_f0(jdim);
            }

        }
        M_I4 = M_fof.contract(M_Ck);
        M_I5 = M_fof.contract(M_Ck * M_Ck);
        M_Cfof = M_Ck * M_fof;
    }
}



void
TransverselyIsoytopicMaterial::evaluateStress( ElasticSolverType solverType)
{
    IsotropicMaterial::evaluateStress(solverType);
}

/// This method is used only for mixed elasticy
void
TransverselyIsoytopicMaterial::evaluateDeviatoricJacobian(  const libMesh::TensorValue <double>&  dU, double q )
{
    IsotropicMaterial::evaluateDeviatoricJacobian(dU, q);
    auto  dF = dU;
    auto  dC = M_Fk.transpose() * dF + dF.transpose() * M_Fk;
    auto CinvdCCinv = M_Cinvk * dC * M_Cinvk;
    auto Jac = 0.0 * M_identity;

    auto dI4 = dC.contract(M_fof);
    auto dI5 = 2.0 * dC.contract(M_Cfof);
    auto dW4 = M_W44 * dI4 + M_W45 * dI5;
    auto dW5 = M_W54 * dI4 + M_W55 * dI5;

    // I4 part
//    //    Jac += 2.0 * dW4 * ( M_fof - M_I4 / 3.0 * M_Cinvk );
//    auto f = M_Fk * M_f0;
//    auto df = dU * M_f0;
//    Jac += 2.0 * dI4 * ( M_fof );
//
//     M_deviatoric_jacobian += M_Fk * Jac;
    Jac += 2.0 * dW4 * ( M_fof - M_I4 / 3.0 * M_Cinvk );
    Jac += 2.0 * M_W4 * ( CinvdCCinv * M_I4 / 3.0 - dI4 / 3.0 * M_Cinvk );

    // I5 part
    auto dCfof = dC * M_fof;

    if(M_active)
    {
        // dCfof = CAinv * dC * fof
        dCfof = M_CAinv * dCfof;
    }

    Jac += 2.0 * dW5 * ( M_Cfof + M_Cfof.transpose() - 2.0 * M_I5 / 3.0 * M_Cinvk);
    Jac += 2.0 * M_W5 * ( dCfof + dCfof.transpose() + CinvdCCinv * M_I5 * 2.0 / 3.0 - 2.0 * dI5 / 3.0 * M_Cinvk );

    M_deviatoric_jacobian += M_Fk * Jac;
}


/// This method is used for primal elasticity
void
TransverselyIsoytopicMaterial::evaluateJacobian(  const libMesh::TensorValue <double>&  dU, double q)
{
    // NOT IMPLEMENTED
    // ONLY INCOMPRESSIBLE
       std::cout << "--- TransverselyIsoytopicMaterial : evaluateJacobian not coded for TransverselyIsoytopicMaterial !!!" << std::endl;
       throw std::runtime_error("evaluateJacobian not coded for  TransverselyIsoytopicMaterial");
}


double
TransverselyIsoytopicMaterial::W4 (double I4, double I5, double J)
{
    double W4 = 0.0;
//    if(I4 < 1) return 0.0;
    switch(M_hexagonalCL)
    {
        case HexagonalCL::Murphy:
        {
            double mu4 = M_parameters[3];
            double mu5 = M_parameters[4];
            W4 = mu5 + mu4 * (I4 - 1);
            break;
        }
        case HexagonalCL::Exponential:
        {
            double af = M_parameters[3];
            double bf = M_parameters[4];
            W4 = af * ( I4 - 1 ) * std::exp( bf * ( I4 - 1 ) * ( I4 - 1 ) );
            break;
        }
        default:
        {
            break;
        }
    }
    return W4;
}

double
TransverselyIsoytopicMaterial::W44(double I4, double I5, double J)
{
    double W44 = 0.0;
//    if(I4 < 1) return 0.0;
    switch(M_hexagonalCL)
    {
        case HexagonalCL::Murphy:
        {
            W44 = M_parameters[3]; // mu4
            break;
        }
        case HexagonalCL::Exponential:
        {
            double af = M_parameters[3];
            double bf = M_parameters[4];
            double aux = bf * ( I4 - 1 ) * ( I4 - 1 );
            W44 = (1 + 2 *  aux ) * af * std::exp( aux );
            break;
        }
        default:
        {
            break;
        }
    }
    return W44;
}

double
TransverselyIsoytopicMaterial::W45(double I4, double I5, double J)
{
    double W45 = 0.0;
//    if(I4 < 1) return 0.0;
    switch(M_hexagonalCL)
    {
        case HexagonalCL::Murphy:
        default:
        {
            break;
        }
    }
    return W45;

}

double
TransverselyIsoytopicMaterial::W5 (double I4, double I5, double J)
{
    double W5 = 0.0;
//    if(I4 < 1) return 0.0;
    switch(M_hexagonalCL)
    {
        case HexagonalCL::Murphy:
        {
            W5 = -0.5 * M_parameters[4];
            break;
        }
        default:
        {
            break;
        }
    }
    return W5;
}

double
TransverselyIsoytopicMaterial::W54(double I4, double I5, double J)
{
    double W54 = 0.0;
//    if(I4 < 1) return 0.0;
    switch(M_hexagonalCL)
    {
        case HexagonalCL::Murphy:
        default:
        {
            break;
        }
    }
    return W54;
}

double
TransverselyIsoytopicMaterial::W55(double I4, double I5, double J)
{
    double W55 = 0.0;
//    if(I4 < 1) return 0.0;
    switch(M_hexagonalCL)
    {
        case HexagonalCL::Murphy:
        default:
        {
            break;
        }
    }
    return W55;
}

} /* namespace BeatIt */
