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
 * Elasticity.hpp
 *
 *  Created on: Oct 19, 2016
 *      Author: srossi
 */

#ifndef SRC_ELASTICITY_ELASTICITY_HPP_
#define SRC_ELASTICITY_ELASTICITY_HPP_

// Include file that defines (possibly multiple) systems of equations.
#include "libmesh/equation_systems.h"
#include "libmesh/linear_solver.h"

#include <memory>
#include "libmesh/getpot.h"
#include "BoundaryConditions/BCHandler.hpp"
#include "Util/SpiritFunction.hpp"
#include "Elasticity/ElasticSolverType.hpp"

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
class MeshRefinement;
template<class T> class DenseMatrix;
template<class T> class DenseVector;
class QGauss;
class MeshBase;
}

namespace BeatIt
{
class BCHandler;
class SpiritFunction;
class Material;

class Elasticity
{
public:
    typedef libMesh::ExodusII_IO EXOExporter;
    typedef libMesh::VTKIO VTKExporter;
    typedef libMesh::ExplicitSystem ParameterSystem;
    typedef libMesh::GMVIO Exporter;

    Elasticity(libMesh::EquationSystems& es, std::string system_name);
    virtual void setup(const GetPot& data, std::string section = "elasticity");
    virtual void setupSystem(std::string section = "elasticity");
    virtual void setupParameters(std::string section = "elasticity");
    virtual void save_exo(
            const std::string& output_filename = "elasticity.exo",
            int step = 0,
            double time = 1.0);
    void write_equation_system(const std::string& es = "elasticity.dat");
    void read_equation_system(const std::string& es = "elasticity.dat");

    void deleteSystems();
    virtual void assemble_residual(
            double dt = 0.0,
            libMesh::NumericVector<libMesh::Number>* activation_ptr = nullptr);
    virtual void assemble_jacobian()
    {
    }
    void newton(
            double dt = 0.0,
            libMesh::NumericVector<libMesh::Number>* activation_ptr = nullptr);
    virtual void update_displacements(double /* dt */)
    {
    }
    void solve_system();
    virtual void project_pressure();

    virtual void advance()
    {
    }
    virtual void init_exo_output(const std::string& output_filename);
    virtual void save(const std::string& output_filename, int step);
    virtual ~Elasticity();

    struct NewtonData
    {
        NewtonData()
                : tol(1e-9), max_iter(20), iter(0)
        {
        }
        double tol;
        int max_iter;
        int iter;
    };

    void setTime(double time);

    GetPot M_datafile;
    libMesh::EquationSystems& M_equationSystems;
    BCHandler M_bch;

    std::set<std::string> M_parametersExporterNames;
    std::unique_ptr<EXOExporter> M_exporter;
    std::unique_ptr<Exporter> M_GMVexporter;
    std::unique_ptr<VTKExporter> M_VTKexporter;
    libMesh::UniquePtr<libMesh::LinearSolver<libMesh::Number> > M_linearSolver;
    libMesh::UniquePtr<libMesh::LinearSolver<libMesh::Number> > M_projectionsLinearSolver;
    std::string M_outputFolder;
    SpiritFunction M_rhsFunction;
    std::string M_myName;

    typedef std::unique_ptr<Material> MaterialPtr;
    std::map<unsigned int, MaterialPtr> M_materialMap;

    ElasticSolverType M_solverType;
    NewtonData M_newtonData;
    bool M_stabilize;

    void evaluate_nodal_I4f();

protected:
    void apply_BC(
            const libMesh::Elem*& elem,
            libMesh::DenseMatrix<libMesh::Number>& Ke,
            libMesh::DenseVector<libMesh::Number>& Fe,
            libMesh::UniquePtr<libMesh::FEBase>& fe_face,
            libMesh::QGauss& qface,
            const libMesh::MeshBase& mesh,
            int n_ux_dofs,
            MaterialPtr mat = nullptr,
            double dt = 0.0);
};

} /* namespace BeatIt */

#endif /* SRC_ELASTICITY_ELASTICITY_HPP_ */
