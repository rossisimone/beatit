// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// <h1>Systems Example 1 - Stokes Equations</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This example shows how a simple, linear system of equations
// can be solved in parallel.  The system of equations are the familiar
// Stokes equations for low-speed incompressible fluid flow.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>

// Basic include file needed for the mesh functionality.
#include "libmesh/mesh_smoother_laplace.h"
#include "libmesh/libmesh.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/boundary_info.h"
#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/explicit_system.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/analytic_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/transient_system.h"

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#include <libmesh/point_locator_tree.h>
#include "libmesh/getpot.h"
#include "Util/IO/io.hpp"

static bool fix_M72 = false;
static std::set<unsigned int> SPV_ids;
static std::set<unsigned int> IPV_ids;
// Bring in everything from the libMesh namespace
using namespace libMesh;




void evaluate_average_mv_translation_and_normal(MeshBase& mesh,
                                                double& nx, double& ny, double& nz,
                                                double& tx, double& ty, double& tz, unsigned int mv_boundary_id = 4)
{

   EquationSystems es(mesh);
   ExplicitSystem& system = es.add_system<ExplicitSystem>("integrate");
   system.add_variable("dummy", FIRST);
   es.init();
   int dim = mesh.mesh_dimension();
   const DofMap & dof_map = system.get_dof_map();
   FEType fe_type = dof_map.variable_type(0);
   std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
   QGauss qface(dim-1, FIRST);
   fe_face->attach_quadrature_rule (&qface);
   const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
   const std::vector<Real> & JxW_face = fe_face->get_JxW();
   const std::vector<Point> & qface_point = fe_face->get_xyz();
   const std::vector<Point> & normals = fe_face->get_normals();

   // evaluate Area 
   double area = 0.0;
   double Nx = 0.0, Ny = 0.0, Nz = 0.0;
   double Tx = 0.0, Ty = 0.0, Tz = 0.0;
   for ( auto & elem : mesh.element_ptr_range() )
   {
        for (int side = 0; side < elem->n_sides(); side++)
        {
            boundary_id_type bid = mesh.boundary_info->boundary_id(elem, side);
            if(bid >= 4 && fix_M72)
            {
              mesh.boundary_info->remove_side(elem, side, bid);
              mesh.boundary_info->add_side(elem, side, bid+1);
            }
            if( elem->neighbor_ptr(side) == nullptr && fix_M72 )
            {
              bool assign_mv = true;
              for (int k = 2; k <= 11; ++k) 
              {
                 if(k == bid)
                 {
                     assign_mv = false;
                     break;
                  }
               }
              if(assign_mv) mesh.boundary_info->add_side(elem, side, 4);
            }
              
            if ( mesh.boundary_info->has_boundary_id( elem, side, mv_boundary_id))
            {
              fe_face->reinit(elem, side);       
              for (int qp = 0; qp <= qface.n_points(); qp++)
              {
                  area += JxW_face[qp]; 
                  Tx += JxW_face[qp] * qface_point[qp](0);
                  Ty += JxW_face[qp] * qface_point[qp](1);
                  Tz += JxW_face[qp] * qface_point[qp](2);
                  Nx += JxW_face[qp] * normals[qp](0);
                  Ny += JxW_face[qp] * normals[qp](1);
                  Nz += JxW_face[qp] * normals[qp](2);
              }
            }
        }
   }
   std::cout << "Area: " << area << std::endl;
   nx = Nx / area;
   ny = Ny / area;
   nz = Nz / area;
   tx = Tx / area;
   ty = Ty / area;
   tz = Tz / area;
}

void evaluate_pv_direction(MeshBase& mesh, double& tx, double& ty, double& tz, std::set<unsigned int> SPV_ids, std::set<unsigned int> IPV_ids)
{

   EquationSystems es(mesh);
   ExplicitSystem& system = es.add_system<ExplicitSystem>("integrate");
   system.add_variable("dummy", FIRST);
   es.init();
   int dim = mesh.mesh_dimension();
   const DofMap & dof_map = system.get_dof_map();
   FEType fe_type = dof_map.variable_type(0);
   std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
   QGauss qface(dim-1, FIRST);
   fe_face->attach_quadrature_rule (&qface);
   const std::vector<std::vector<Real>> & phi_face = fe_face->get_phi();
   const std::vector<Real> & JxW_face = fe_face->get_JxW();
   const std::vector<Point> & qface_point = fe_face->get_xyz();
   const std::vector<Point> & normals = fe_face->get_normals();
   // evaluate Area
   double areaIPV = 0.0;
   double areaSPV = 0.0;
   double Nx = 0.0, Ny = 0.0, Nz = 0.0;
   double P1x = 0.0,P1y = 0.0,P1z = 0.0;
   double P2x = 0.0,P2y = 0.0,P2z = 0.0;
   for ( auto & elem : mesh.element_ptr_range() )
   {
        for (int side = 0; side < elem->n_sides(); side++)
        {
            bool is_on_SPV = false;
            bool is_on_IPV = false;
            for(auto & id : SPV_ids)
            {
               is_on_SPV = mesh.boundary_info->has_boundary_id( elem, side, id);
               if(is_on_SPV) break;
            }
            if(!is_on_SPV)
            {
                for(auto & id : IPV_ids)
                {
                   is_on_IPV = mesh.boundary_info->has_boundary_id( elem, side, id);
                   if(is_on_IPV) break;
                }
            }
            if ( is_on_IPV || is_on_SPV)
            {
              fe_face->reinit(elem, side);
              for (int qp = 0; qp <= qface.n_points(); qp++)
              {
                  if(is_on_IPV)
                  {
                     areaIPV += JxW_face[qp];
                     P1x += JxW_face[qp] * qface_point[qp](0);
                     P1y += JxW_face[qp] * qface_point[qp](1);
                     P1z += JxW_face[qp] * qface_point[qp](2);
                  }
                  else
                  {
                     areaSPV += JxW_face[qp];
                     P2x += JxW_face[qp] * qface_point[qp](0);
                     P2y += JxW_face[qp] * qface_point[qp](1);
                     P2z += JxW_face[qp] * qface_point[qp](2);
                  }
              }
            }
        }
   }
   P1x /= areaIPV;
   P1y /= areaIPV;
   P1z /= areaIPV;
   P2x /= areaSPV;
   P2y /= areaSPV;
   P2z /= areaSPV;

   tx = P2x - P1x;
   ty = P2y - P1y;
   tz = P2z - P1z;
   std::cout << tx << " " << ty << " " << tz<< std::endl;
}


void rotate_mesh_to_xy_plane(MeshBase& mesh, double scale = 1.0)
{
    double nx, ny, nz, tx, ty, tz;
    evaluate_average_mv_translation_and_normal(mesh, nx, ny, nz, tx, ty, tz);
    std::cout << "translating: " << tx << ", " << ty << ", " << tz << std::endl;
    std::cout << "rotating: " << nx << ", " << ny << ", " << nz << std::endl;
    std::cout << "scaling: " << scale << std::endl;
    libMesh::Point n(nx,ny, nz);
    n = n.unit();
    libMesh::Point n_xy(0, 0, 1.0);
    libMesh::Point n_x(1, 0, 0.0);
    // rotation axis
    std::cout << "Rotation axis: ";
    libMesh::Point u = n_xy.cross(n);
    u = u.unit();
    u.print(std::cout);
    std::cout << std::endl;
    // rotation angle
    double cos_tetha = n * n_xy;
    double tetha = -std::acos(cos_tetha);
    std::cout << "Rotation angle " << tetha << std::endl; 
    double sin_tetha = std::sin(tetha);
    // translation 
    libMesh::Point t(tx, ty,tz);

    TensorValue<Number> R;
    R(0,0) =          cos_tetha + u(0) * u(0) * ( 1 - cos_tetha);
    R(0,1) = - u(2) * sin_tetha + u(0) * u(1) * ( 1 - cos_tetha);
    R(0,2) =   u(1) * sin_tetha + u(0) * u(2) * ( 1 - cos_tetha);
    R(1,0) =   u(2) * sin_tetha + u(1) * u(0) * ( 1 - cos_tetha);
    R(1,1) =          cos_tetha + u(1) * u(1) * ( 1 - cos_tetha);
    R(1,2) = - u(0) * sin_tetha + u(1) * u(2) * ( 1 - cos_tetha);
    R(2,0) = - u(1) * sin_tetha + u(2) * u(0) * ( 1 - cos_tetha);
    R(2,1) =   u(0) * sin_tetha + u(2) * u(1) * ( 1 - cos_tetha);
    R(2,2) =          cos_tetha + u(2) * u(2) * ( 1 - cos_tetha);
/*    for(auto & node : mesh.node_ptr_range() )
    {
        const libMesh::Point X = *node;
        libMesh::Point Tx = X - t;
        libMesh::Point STx = scale * Tx;
        libMesh::Point uxSTx = u.cross(STx);
        *node = Point( RSTx(0), RSTx(1), RSTx(2));
    }
*/
    for(auto & node : mesh.node_ptr_range() )
    {
        const libMesh::Point X = *node;
        libMesh::Point Tx = X - t;
        libMesh::Point STx = scale * Tx;
        libMesh::Point uxSTx = u.cross(STx);
        //libMesh::Point RSTx =  u * ( u * STx ) + cos_tetha * uxSTx.cross(u) + sin_tetha * uxSTx;
        libMesh::Point RSTx =  R * STx;
        // apply rotation 180 wrt x-axis
        *node = Point( RSTx(0), -RSTx(1), -RSTx(2));
    }

    evaluate_pv_direction(mesh, tx, ty, tz, SPV_ids, IPV_ids);
    Point pv_direction(tx, ty, 0);
    pv_direction = pv_direction.unit();
    pv_direction.print(std::cout); 
    u = n_xy;
    cos_tetha = pv_direction * n_x;
    tetha = -std::acos(cos_tetha);
    std::cout << "\nX Rotation angle " << tetha << std::endl;
    sin_tetha = std::sin(tetha);
    // translation
    R(0,0) =          cos_tetha + u(0) * u(0) * ( 1 - cos_tetha);
    R(0,1) = - u(2) * sin_tetha + u(0) * u(1) * ( 1 - cos_tetha);
    R(0,2) =   u(1) * sin_tetha + u(0) * u(2) * ( 1 - cos_tetha);
    R(1,0) =   u(2) * sin_tetha + u(1) * u(0) * ( 1 - cos_tetha);
    R(1,1) =          cos_tetha + u(1) * u(1) * ( 1 - cos_tetha);
    R(1,2) = - u(0) * sin_tetha + u(1) * u(2) * ( 1 - cos_tetha);
    R(2,0) = - u(1) * sin_tetha + u(2) * u(0) * ( 1 - cos_tetha);
    R(2,1) =   u(0) * sin_tetha + u(2) * u(1) * ( 1 - cos_tetha);
    R(2,2) =          cos_tetha + u(2) * u(2) * ( 1 - cos_tetha);

    for(auto & node : mesh.node_ptr_range() )
    {
        const libMesh::Point RSTx = *node;
        libMesh::Point RRSTx =  R * RSTx;
        *node = Point( RRSTx(0), RRSTx(1), RRSTx(2));
    }
}
// The main program.
int main(int argc, char ** argv)
{
    // Initialize libMesh.
    LibMeshInit init(argc, argv);
    // Read input file
    GetPot commandLine ( argc, argv );
    std::string datafile_name = commandLine.follow ( "data.beat", 2, "-i", "--input" );
    GetPot data(datafile_name);

    std::string mesh_file = data("mesh", "NONE");
    ParallelMesh mesh(init.comm());
    mesh.read(mesh_file + ".msh");
    mesh.print_info();
    if(mesh_file == "M72") fix_M72 = true;
    bool export_exodus = data("export_exodus", true);
    bool export_vtk = data("export_vtk", true);
    bool export_nemesis = data("export_nemesis", true);
    double scale = data("scale", 1.0);
    //if(nz != 1.0)
    {
          //evaluate_average_mv_translation_and_normal(mesh, nx, ny, nz, tz, ty, tz);
          std::string  SPV_ids_list = data("SPV_ids", "");
          std::string  IPV_ids_list = data("IPV_ids", "");
          BeatIt::readList(SPV_ids_list, SPV_ids);
          BeatIt::readList(IPV_ids_list, IPV_ids);
          std::cout << "Rotating ... " << std::endl;
          rotate_mesh_to_xy_plane(mesh, scale);
    }
    //else
    //{
    //      std::cout << "nz: " << nz << std::endl;
    //}

    // evaluate mesh statistics
    bool print_mesh_statistics = data("info",true);
    if(print_mesh_statistics)
    {
    std::vector<std::pair< dof_id_type, dof_id_type> > node_map;
    for (auto & elem : mesh.element_ptr_range() )
    {
        for(unsigned int edge = 0; edge < elem->n_edges(); ++edge)
        {
             std::unique_ptr<libMesh::Elem> edge_ptr = elem->build_edge_ptr(edge);
             auto p1id = edge_ptr->node_ptr(0)->id();
             auto p2id = edge_ptr->node_ptr(1)->id();
             std::pair<dof_id_type, dof_id_type> node_id_pair(p1id, p2id);
             auto p1 = std::find(node_map.begin(), node_map.end(), node_id_pair);
             std::pair<dof_id_type, dof_id_type> node_id_pair_reversed(p2id, p1id);
             auto p2 = std::find(node_map.begin(), node_map.end(), node_id_pair_reversed);
             if( p1 == node_map.end() && p2 == node_map.end())
                       node_map.push_back(node_id_pair);
        }
    }
    int nedges = node_map.size();
    std::cout << "create node map with size: " << nedges << std::endl;

    std::vector<double> hvec(nedges);
    double h_min = 1000000.0, h_max = 0, h_average = 0, h_median = 0.0;
    for(int k = 0; k < nedges; ++k )
    {
        Point p1(mesh.node_ref(node_map[k].first));
        Point p2(mesh.node_ref(node_map[k].second));
        p2 -= p1;
        hvec[k] = p2.norm();

        h_min = std::min(h_min, hvec[k]);
        h_max = std::max(h_max, hvec[k]);
        h_average += hvec[k];
    }
    h_average /= nedges;

    {
        size_t n = hvec.size() / 2;
        std::nth_element(hvec.begin(), hvec.begin()+n, hvec.end());
        h_median = hvec[n];
    }
    for(int nref = 0; nref < 5; nref++)
    {
        double div = std::pow(2,nref);
        std::cout << "Mesh " << nref << ": h_min: " << h_min/div << ", h_max: " << h_max/div  << ", h_mean: " << h_average/div  << ", h_median: " << h_median/div << std::endl;
    }
    }
    // Refine mesh
    std::string output_name = data("output", mesh_file);
    if(export_exodus) ExodusII_IO (mesh).write(output_name+"_m0.e");
    if(export_vtk) VTKIO (mesh).write(output_name+"_m0.pvtu");
    if(export_nemesis) Nemesis_IO (mesh).write(output_name+"_m0.nem");

    int nref = data("nref", 0);
    int nsm = data("nsm", 0);
    if(nref > 0)
    {
      for(int nr = 1; nr <= nref; nr++)
      {
          ParallelMesh refined_mesh(mesh);
          MeshRefinement refine(refined_mesh);
          std::cout << "Refine: " << nr << std::endl;
          refine.uniformly_refine(nr);
          //mesh.prepare_for_use(true);
          std::cout << "Exporting: " << nr << std::endl;
          if(export_exodus) ExodusII_IO (refined_mesh).write(output_name+"_m"+std::to_string(nr)+".e");
          if(export_vtk) VTKIO (refined_mesh).write(output_name+"_m"+std::to_string(nr)+".pvtu");
          if(export_nemesis) Nemesis_IO (refined_mesh).write(output_name+"_m"+std::to_string(nr)+".nem");
      //if(nsm)
      //{
      //    LaplaceMeshSmoother smoother(mesh);
      //    smoother.init();
      //    std::cout << "smoothing " << nsm << " times " << std::endl; 
      //    smoother.smooth(nsm);
      //}
      }
    }

    return 0;
}
