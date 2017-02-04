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
 * \file Electromechanics.hpp
 *
 * \class Electromechanics
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
 * Created on: Dec 23, 2016
 *
 */

#ifndef SRC_ELECTROMECHANICS_ELECTROMECHANICS_HPP_
#define SRC_ELECTROMECHANICS_ELECTROMECHANICS_HPP_

#include "Elasticity/MixedElasticity.hpp"
#include "Electrophysiology/Monodomain/Monowave.hpp"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"


#include <string>
#include <memory>

// Forward Declarations
namespace libMesh
{
class Mesh;
class ExodusII_IO;
class VTKIO;
class GMVIO;
class MeshRefinement;
class TimeData;
class ExplicitSystem;
class LinearImplicitSystem;
class PacingProtocol;
class ErrorVector;
class  MeshRefinement;
class MeshBase;
}



namespace BeatIt
{

class Elasticity;
class MixedElasticity;
class Monowave;
/*!
 *
 */
class Electromechanics
{
public:
    Electromechanics( libMesh::EquationSystems& es,
                      std::string system_name );
    ~Electromechanics( );
    /// The datafile will contain only the addresses to the specific seprate datafiles
    void setup( GetPot& data,
                std::string electro_section = "monodomain",
                std::string elasticity_section = "elasticity",
                std::string activation_section = "activation");

    void init(double time);

    void compute_activation(double dt);

    GetPot                      M_datafile;
    std::unique_ptr<MixedElasticity> M_elasticity;
    std::unique_ptr<Monowave>   M_monowave;
    libMesh::EquationSystems&  M_equationSystems;
    std::string M_myName;
    std::string M_outputFolder;
    libMesh::MeshBase& M_mesh;

    typedef libMesh::ExodusII_IO EXOExporter;
    std::unique_ptr<EXOExporter> M_exporter;
    void save_exo(int step, double time);

    void solve_mechanics();

    void solve_reaction_step( double dt,
                              double time,
                              int step = 0,
                              bool useMidpoint = true,
                              const std::string& mass = "mass");

    void assemble_electrophysiology_matrices();

};

} /* namespace BeatIt */

#endif /* SRC_ELECTROMECHANICS_ELECTROMECHANICS_HPP_ */
