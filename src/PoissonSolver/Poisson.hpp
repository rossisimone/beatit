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
 * Poisson.hpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */

#ifndef SRC_POISSONSOLVER_POISSON_HPP_
#define SRC_POISSONSOLVER_POISSON_HPP_

// Include file that defines (possibly multiple) systems of equations.
#include "libmesh/equation_systems.h"
#include "libmesh/linear_solver.h"
#include "libmesh/petsc_linear_solver.h"

#include <memory>
#include "libmesh/getpot.h"
#include "BoundaryConditions/BCHandler.hpp"
#include "Util/SpiritFunction.hpp"


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
template < class T > class DenseMatrix;
template < class T > class DenseVector;
class QGauss;
class MeshBase;
}


namespace BeatIt {

class BCHandler;
class SpiritFunction;

class Poisson {
	    typedef libMesh::ExodusII_IO EXOExporter;

public:
//	Poisson( libMesh::MeshBase & mesh );
    Poisson( libMesh::EquationSystems& es, std::string system_name = "poisson" );
	virtual ~Poisson();
    void setup(const GetPot& data, std::string section = "poisson" );
    void assemble_system();
    void solve_system();
    void save_exo(const std::string& output_filename = "poisson.exo");


    const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
    get_gradient();

    const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
    get_solution();

    const libMesh::UniquePtr<libMesh::NumericVector<libMesh::Number> >&
    get_P0_solution();

    void write_equation_system(const std::string& es = "poisson.dat");
    void read_equation_system(const std::string& es = "poisson.dat");

    double get_solution_norm();

    void compute_elemental_solution_gradient();

    void deleteSystems();
    /// input file

    GetPot                     M_datafile;
    libMesh::EquationSystems&  M_equationSystems;
    BCHandler M_bch;

    std::set<std::string> M_parametersExporterNames;
    std::unique_ptr<EXOExporter> M_exporter;
    std::unique_ptr<libMesh::PetscLinearSolver<libMesh::Number> > M_linearSolver;
//    libMesh::UniquePtr<libMesh::LinearSolver<libMesh::Number> > M_linearSolver;
    std::string  M_outputFolder;
    SpiritFunction M_rhsFunction;
    std::string M_myName;
    std::string M_myNameGradient;
    std::string M_myNameP0;

private:
    void apply_BC( const libMesh::Elem*& elem,
            libMesh::DenseMatrix<libMesh::Number>& Ke,
            libMesh::DenseVector<libMesh::Number>& Fe,
            libMesh::UniquePtr<libMesh::FEBase>& fe_face,
            libMesh::QGauss& qface,
            const libMesh::MeshBase& mesh);

};

} /* namespace BeatIt */

#endif /* SRC_POISSONSOLVER_POISSON_HPP_ */
