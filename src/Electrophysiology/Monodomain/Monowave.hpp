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
 * Monowave.hpp
 *
 *  Created on: Oct 27, 2016
 *      Author: srossi
 */

#ifndef SRC_ELECTROPHYSIOLOGY_MONODOMAIN_MONOWAVE_HPP_
#define SRC_ELECTROPHYSIOLOGY_MONODOMAIN_MONOWAVE_HPP_

#include "Electrophysiology/ElectroSolver.hpp"

namespace BeatIt
{


/// Class
class Monowave : public virtual ElectroSolver
{
public:
    typedef libMesh::GMVIO Exporter;
    // Another alternative when not using AMR
    typedef libMesh::ExodusII_IO EXOExporter;

    /// Empty construcor
//    Monodomain( libMesh::MeshBase & mesh );
    Monowave( libMesh::EquationSystems& es );
    ~Monowave();
    void setup_systems(GetPot& data, std::string section = "monowave" );

    void init_systems(double time);

    void amr( libMesh:: MeshRefinement& mesh_refinement, const std::string& type = "kelly" );


    void cut(double time, std::string f);
    void assemble_matrices(double dt = 1.0);
    void assemble_dg_matrices(double dt = 1.0);
    void assemble_cg_matrices(double dt = 1.0);
    void setup_local_conductivity(libMesh::TensorValue<double>& D0,
                                           double Dff, double Dss, double Dnn,
                                           double * f0, double * s0, double * n0);
    void form_system_matrix(double dt, bool useMidpoint = true, const std::string& mass = "lumped_mass");
    void form_system_rhs(double dt, bool useMidpoint = true, const std::string& mass = "lumped_mass");

//    void solve_reaction_step( double dt,
//                              double time,
//                              int step = 0,
//                              bool useMidpoint = true,
//                              const std::string& mass = "mass",
//                              libMesh::NumericVector<libMesh::Number>* I4f_ptr = nullptr);

    void solve_diffusion_step(double dt, double time,  bool useMidpoint = true, const std::string& mass = "lumped_mass", bool reassemble = true);

    void generate_fibers(   const GetPot& data,
                            const std::string& section = "rule_based_fibers" );

};


ElectroSolver* createMonowave( libMesh::EquationSystems& es );

namespace
{
    static bool register_monowave = BeatIt::ElectroSolver::ElectroFactory::Register("monowave", &createMonowave);
}

} /* namespace BeatIt */


#endif /* SRC_ELECTROPHYSIOLOGY_MONODOMAIN_MONOWAVE_HPP_ */
