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
 * .hpp
 *
 *  Created on: Oct 27, 2016
 *      Author: srossi
 */

#ifndef SRC_ELECTROPHYSIOLOGY_Bidomain_HPP_
#define SRC_ELECTROPHYSIOLOGY_Bidomain_HPP_

#include "Electrophysiology/ElectroSolver.hpp"

namespace BeatIt
{



//enum class BidomainAnisotropy {Isotropic, TransverselyIsotropic, Orthotropic };
//enum class BidomainEquationType { ParabolicEllipticBidomain,
//	                              ParabolicEllipticHyperbolic,
//								  ParabolicParabolicHyperbolic };
//
//enum class BidomainTimeIntegrator { FirstOrderIMEX, SecondOrderIMEX };

/// Class
class BidomainWithBath : public virtual ElectroSolver
{
public:
    /// Empty construcor
//    Monodomain( libMesh::MeshBase & mesh );
    BidomainWithBath( libMesh::EquationSystems& es );
    ~BidomainWithBath();
    void setup_systems(GetPot& data, std::string section = "");

    void init_systems(double time);
    void assemble_matrices(double dt = 1.0);
    void form_system_matrix(double dt, bool useMidpoint = true, const std::string& mass = "lumped_mass") {}
    void form_system_rhs(double dt, bool useMidpoint = true, const std::string& mass = "lumped_mass");
//    void solve_reaction_step( double dt,
//                              double time,
//                              int step = 0,
//                              bool useMidpoint = true,
//                              const std::string& mass = "mass",
//                              libMesh::NumericVector<libMesh::Number>* I4f_ptr = nullptr);

    void solve_diffusion_step(double dt, double time,  bool useMidpoint = true, const std::string& mass = "lumped_mass", bool reassemble = true);



    // For Dirichlet Boundary conditions
    // We use the boundary info to get the list of node with the specified ID
    std::vector< libMesh::dof_id_type >      M_node_id_list;
    std::vector< libMesh::boundary_id_type > M_bc_id_list;
    // For ground node
    int M_ground_point_id;
};



ElectroSolver* createBidomainWithBath( libMesh::EquationSystems& es );

namespace
{
    static bool register_bidomain_bath = BeatIt::ElectroSolver::ElectroFactory::Register("bidomainbath", &createBidomainWithBath);
}


} /* namespace BeatIt */


#endif /* SRC_ELECTROPHYSIOLOGY_Bidomain_HPP_ */
