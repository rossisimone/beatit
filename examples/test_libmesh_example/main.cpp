// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// <h1>Introduction Example 3 - Solving a Poisson Problem</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This is the third example program.  It builds on
// the second example program by showing how to solve a simple
// Poisson system.  This example also introduces the notion
// of customized matrix assembly functions, working with an
// exact solution, and using element iterators.
// We will not comment on things that
// were already explained in the second example.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Use input file
#include "libmesh/getpot.h"

// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/vtk_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/equation_systems.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/enum_solver_package.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Include headerfile to detect boundary id
#include "libmesh/boundary_info.h"

// To read the IDs from the input file
#include "Util/IO/io.hpp"

// Bring in everything from the libMesh namespace
using namespace libMesh;


// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the  EquationSystems object and the
// name of the system we are assembling as input.  From the
//  EquationSystems object we have access to the  Mesh and
// other objects we might need.
void assemble_poisson(EquationSystems& es,
    const std::string& system_name);

// Function prototype for the exact solution.
Real exact_solution(const Real x,
    const Real y,
    const Real z = 0.)
{
    static const Real pi = acos(-1.);

    //  return (16/9)*cos(.25*pi*(3*x-5))*sin(pi*y)*cos(.5*pi*z);
    //  return cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z); //original
    return sin(.5 * pi * x) * cos(.5 * pi * y) * cos(.5 * pi * z); // reversed - test
}

class Poisson{
	public:
		// Member variable Dirichlet BC list - can use inside poisson_assemble
		// each list stores subregion IDs to where that boundary condition should be imposed
		std::vector<int> dirichlet_side_set_list;
		std::vector<int> dirichlet_bc_list;


		// Stores a list of the names of all input files for each problem we are solving
		std::vector<std::string> input_list;

		// Pointer to a vector of points
		// These vectors store coordinates of centers of spheres and their radius
		// example: we want to impose Dirichlet BC 0.5 on spheres w/ centers at Point1=(0,2,3) and radius 0.3
		//          and another sphere w/ center at Point2=(1,9,8) and radius 0.4
		// 			then dirichlet_id05_points = [Point1 Point2];
		//          and  dirichlet_id05_radius = [0.3 0.4]

		std::vector<libMesh::Point> dirichlet_id_points;

		// Radius of sphere (neighborhood) where we impose Dirichlet BC
		std::vector<double> dirichlet_id_radius;


		void reinit(std::string, std::string,std::string,
				 std::string, std::string, std::string);
};

/*
class BCs{
	id of side sets
	values of bcs
	do i want to apply on the side set or sphere? boolean that depends on radius
	point and radius
	variable negative if wna to apply to a set id and positive to a sphere
};*/


	void Poisson::reinit( std::string dirichlet_x_list, std::string dirichlet_y_list, std::string dirichlet_z_list,
						 std::string dirichlet_r_list, std::string dirichlet_side_set_list_aux, std::string dirichlet_bc_list_aux)
	{

        //erase what was stored in member variables
        dirichlet_side_set_list.clear();
        dirichlet_bc_list.clear();

        dirichlet_id_points.clear();

        dirichlet_id_radius.clear();

        //Read the BC list for the current problem
        BeatIt::readList(dirichlet_side_set_list_aux, dirichlet_side_set_list);
        BeatIt::readList(dirichlet_bc_list_aux, dirichlet_bc_list);

        //for general
        std::vector<double>  x, y, z, r;
        BeatIt::readList(dirichlet_x_list ,     x);
        BeatIt::readList(dirichlet_y_list ,     y);
        BeatIt::readList(dirichlet_z_list ,     z);
        BeatIt::readList(dirichlet_r_list ,     r);

        for (int jj = 0; jj < r.size(); jj++) {
            dirichlet_id_points.emplace_back(x[jj], y[jj], z[jj]);
            dirichlet_id_radius.push_back(r[jj]);
        }

        //check for error
        if (x.size()!=y.size() || x.size()!=z.size()|| x.size()!=r.size())
        {
           	std::cout << "the sizes of x, y, z, and r lists should be the same "<< std::endl;
           	std::cout<< "x.size() = " << x.size()  <<std::endl;
           	std::cout<< "y.size() = " << y.size()  <<std::endl;
           	std::cout<< "z.size() = " << z.size()  <<std::endl;
           	std::cout<< "r.size() = " << r.size()  <<std::endl;
           	throw std::runtime_error("the sizes of x, y, z, and r lists should be the same");
        }

        if(dirichlet_bc_list.size() != dirichlet_side_set_list.size() || dirichlet_bc_list.size() != dirichlet_id_radius.size()){
        	std::cout << "Size of Idirichlet_bc_list and Idirichlet_side_set_list have to be equal"<<std::endl;
        	std::cout << "dirichlet_bc_list.size() = " << dirichlet_bc_list.size() << ", dirichlet_side_set_list.size() = " <<
        			dirichlet_side_set_list.size() << " , dirichlet_id_radius.size() = " << dirichlet_id_radius.size()<<std::endl;
        	throw std::runtime_error("Size of Idirichlet_bc_list and Idirichlet_side_set_list has to be equal");
        }
	}

	//global
	Poisson * poisson_ptr = nullptr;



// MAIN
int main(int argc, char** argv)
{

    // Read input file
    GetPot commandLine(argc, argv);
    std::string datafile_name = commandLine.follow("data.beat", 2, "-i", "--input");
    GetPot data(datafile_name);


    // Modifying code to solve all problems at once
    int N= data("number_of_problems", 0);


    //Create an object of class Poisson and a pointer to it
    Poisson poisson_solver;
    poisson_ptr = & poisson_solver;


    // Read name of the input files corresponding to each problem
    std::string input_list_aux = data("inputs", "");
    BeatIt::readList(input_list_aux, poisson_ptr->input_list);


    // Initialize libraries, like in example 2.
    LibMeshInit init(argc, argv);

    // This example requires a linear solver package.
    libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
        "--enable-petsc, --enable-trilinos, or --enable-eigen");

    // Brief message to the user regarding the program name
    // and command line arguments.
    libMesh::out << "Running " << argv[0];

    for (int i = 1; i < argc; i++)
        libMesh::out << " " << argv[i];

    libMesh::out << std::endl << std::endl;

    // Skip this 2D example if libMesh was compiled as 1D-only.
    libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

    // Create a mesh, with dimension to be overridden later, distributed
    // across the default MPI communicator.
    Mesh mesh(init.comm());

    // read from input
    // Get the ids where Dirichlet BC are imposed
    std::string ids_dirichlet_test = data("dirichletbc", "");  //> did not work
    std::cout << ids_dirichlet_test << " ";
    //ids_dirichlet_bcs=ids_dirichlet_bcs*2;
    int size_dirichlet_bcs = data.vector_variable_size("dirichletbc");
    std::vector<int> ids_dirichlet_bcs; //(size_dirichlet_bcs);
    for (int i = 0; i < size_dirichlet_bcs; i++) {
        //ids_dirichlet_bcs[i]=data("dirichletbc", i, -1);
        //std::cout << data("dirichletbc", i, -1)<< " ";
        ids_dirichlet_bcs.push_back(data("dirichletbc", i, -1));
    }

    // BeatIt::readList(ids_dirichlet_test, ids_dirichlet_bcs);               // IMPORTANT LINE HERE
     //std::cout <<"ids dirichlet:"<< ids_dirichlet_bcs << "\n";
    for (int i = 0; i < size_dirichlet_bcs; i++) {
        std::cout << ids_dirichlet_bcs[i] << " ";
    }

    // Read problem solution index and thresholds
    std::string threshold_pv1 = data("pv1", " 1, 1, 0");
    std::string threshold_pv2 = data("pv2", "1, 1, 0 ");
    std::string threshold_pv3 = data("pv3", "1, 1, 0 ");
    std::string threshold_pv4 = data("pv4", "1, 1, 0");
    std::string threshold_floor = data("floor", "1, 1, 0 ");
    std::string threshold_laa = data("laa", "1, 1, 0");
    std::string threshold_antra1 = data("antra1", "1, 1, 0");
    std::string threshold_antra2 = data("antra2", "1, 1, 0");
    std::string threshold_septum = data("septum", "1, 1, 0");
    std::string threshold_anterior = data("anterior", "1, 1, 0 ");
    std::string threshold_posterior = data("posterior", "1, 1, 0 ");
    std::string threshold_carina1 = data("carina1", "1, 1, 0 ");
    std::string threshold_carina2 = data("carina2", "1, 1, 0 ");

    std::vector<double> tpv1, tpv2, tpv3, tpv4, tfloor, tlaa, tantra1, tantra2, tseptum, tanterior, tposterior, tcarina1, tcarina2;
    BeatIt::readList(threshold_pv1, tpv1);
    BeatIt::readList(threshold_pv2, tpv2);
    BeatIt::readList(threshold_pv3, tpv3);
    BeatIt::readList(threshold_pv4, tpv4);
    BeatIt::readList(threshold_floor,tfloor );
    BeatIt::readList(threshold_laa,tlaa );
    BeatIt::readList(threshold_antra1, tantra1);
    BeatIt::readList(threshold_antra2, tantra2);
    BeatIt::readList(threshold_septum, tseptum);
    BeatIt::readList(threshold_anterior, tanterior);
    BeatIt::readList(threshold_posterior, tposterior);
    BeatIt::readList(threshold_carina1, tcarina1);
    BeatIt::readList(threshold_carina2, tcarina2);

    //Mesh as input?
    bool ifmesh = data("mesh_input", false);

    if (ifmesh == true) {
        // read mesh
        std::string name_mesh = data("mesh", "");
        //mesh.read("F65_m1.e"); // give the full path here if not in the same folder
        mesh.read(name_mesh); // give the full path here if not in the same folder
    }
    else {
        // Use the MeshTools::Generation mesh generator to create a uniform
        // 2D grid on the square [-1,1]^2.  We instruct the mesh generator
        // to build a mesh of 15x15 QUAD9 elements.  Building QUAD9
        // elements instead of the default QUAD4's we used in example 2
        // allow us to use higher-order approximation.
        MeshTools::Generation::build_square(mesh,
            15, 15,
            -1., 1.,
            -1., 1.,
            QUAD9);
    }
    // Print information about the mesh to the screen.
    // Note that 5x5 QUAD9 elements actually has 11x11 nodes,
    // so this mesh is significantly larger than the one in example 2.
    mesh.print_info();


    // Create an equation systems object.
    EquationSystems equation_systems(mesh);

    // Define fiber systems:
    //EquationSystems es(mesh);
    typedef libMesh::ExplicitSystem  FiberSystem;
    FiberSystem& fiber_system = equation_systems.add_system<FiberSystem>("fibers");
    fiber_system.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.init();
    auto& f_v = fiber_system.solution;
    FiberSystem& sheet_system = equation_systems.add_system<FiberSystem>("sheets");
    sheet_system.add_variable("sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
    sheet_system.add_variable("sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
    sheet_system.add_variable("sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
    sheet_system.init();
    auto& s_v = sheet_system.solution;
    FiberSystem& xfiber_system = equation_systems.add_system<FiberSystem>("xfibers");
    xfiber_system.add_variable("xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    xfiber_system.add_variable("xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    xfiber_system.add_variable("xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    xfiber_system.init();
    auto& n_v = xfiber_system.solution;

    // Declare the Poisson system and its variables.
    // The Poisson system is another example of a steady system.
    LinearImplicitSystem& system = equation_systems.add_system<LinearImplicitSystem>("Poisson");

    // Adds the variable "u" to "Poisson".  "u"
    // will be approximated using second-order approximation.
    if (ifmesh == true) {
        equation_systems.get_system("Poisson").add_variable("u", FIRST); //TET4 ELEMENTS
    }
    else {
        equation_systems.get_system("Poisson").add_variable("u", SECOND); // QUAD9 ELEMENTS FOR SQUARE
    }


    // Give the system a pointer to the matrix assembly
    // function.  This will be called when needed by the
    // library.
    equation_systems.get_system("Poisson").attach_assemble_function(assemble_poisson);


    // Initialize the data structures for the equation system.
    equation_systems.init();

    // Prints information about the system to the screen.
    equation_systems.print_info();

    // Reading the input for each of the problems
    for (int i = 0; i < N; i++) {
    	std::cout << "----------------------------------------\n";
    	std::cout << "----------------------------------------\n";

        std::string ui = "u" + std::to_string(i);    //how cool! to use i in the name here!!!
        ExplicitSystem& si = equation_systems.add_system<ExplicitSystem>(ui);
        si.add_variable(ui, FIRST);
        si.init();

        //name of the problem
        GetPot data_i(poisson_ptr->input_list[i]);

        // Read Dirichlet BC ID list from input file
        std::string dirichlet_id_list_aux0 = data_i("dirichletbc0", "");
        std::string dirichlet_id_list_aux1 = data_i("dirichletbc1", "");
        std::string dirichlet_id_list_aux05 = data_i("dirichletbc05", "");
        std::string dirichlet_side_set_list_aux = data_i("Idirichlet_side_set_list", "");
        std::string dirichlet_bc_list_aux = data_i("Idirichlet_bc_list", "");

        // Read Dirichlet BC ID list for septum and laa points
        std::string dirichlet_x0_list = data_i("x0", " ");
        std::string dirichlet_y0_list = data_i("y0", " ");
        std::string dirichlet_z0_list = data_i("z0", " ");
        std::string dirichlet_r0_list = data_i("r0", " ");
        std::string layer0_aux = data_i("Ilayer0", " ");

        std::string dirichlet_x1_list = data_i("x1", " ");
        std::string dirichlet_y1_list = data_i("y1", " ");
        std::string dirichlet_z1_list = data_i("z1", " ");
        std::string dirichlet_r1_list = data_i("r1", " ");
        std::string layer1_aux = data_i("Ilayer1", " ");

 	   std::string dirichlet_x05_list = data_i("x05", " ");
 	   std::string dirichlet_y05_list = data_i("y05", " ");
 	   std::string dirichlet_z05_list = data_i("z05", " ");
 	   std::string dirichlet_r05_list = data_i("r05", " ");
 	   std::string layer05_aux = data_i("Ilayer05", " ");


       std::string dirichlet_x_list = data_i("Ix", " ");
       std::string dirichlet_y_list = data_i("Iy", " ");
       std::string dirichlet_z_list = data_i("Iz", " ");
       std::string dirichlet_r_list = data_i("Ir", " ");

       std::cout << "Ix, Iy, Iz and Ir are respectively: "<< dirichlet_x_list
    		    << " | "<< dirichlet_y_list <<" | "<< dirichlet_z_list <<
				" | "<< dirichlet_r_list << std::endl;


       // reinit Poisson member variables using input file values
 	   poisson_solver.reinit(  dirichlet_x_list,  dirichlet_y_list,  dirichlet_z_list,
							   dirichlet_r_list,  dirichlet_side_set_list_aux, dirichlet_bc_list_aux);

       //Check
       if(poisson_ptr->dirichlet_id_points.size()!= poisson_ptr->dirichlet_id_radius.size()){
               	std::cout << " The number of points does not match the number of radii" << std::endl;
               	std::cout << "for problem i= " << i << std::endl;
               	std::cout << "the size of vector points is "<< poisson_ptr->dirichlet_id_points.size()<<std::endl;
               	std::cout << "the size of vector radius is "<< poisson_ptr->dirichlet_id_radius.size()<<std::endl;
               	throw std::runtime_error(" The number of points does not match the number of radii");
       }


	   //Print info on screen
        std::cout << "Setting Dirichlet BCs on " << poisson_ptr->dirichlet_id_points.size() << " points" << std::endl;

        /*
        // CHANGED LINE 288 - not attaching automatically
        equation_systems.get_system("Poisson").matrix->zero();
        equation_systems.get_system("Poisson").rhs->zero();
        //Assemble function manually
        assemble_poisson( ); // now I can change the inputs here!

        //  Libmesh zeros out matrix and RHS when i call solve
        // preventing that to happen:
        equation_systems.get_system("Poisson").zero_out_matrix_and_rhs=false;*/

        equation_systems.get_system("Poisson").solve();

        // Copy the content of "u" (solution to Poisson problem) to "ui"
        *si.solution = *equation_systems.get_system("Poisson").solution; // "unique vector"
        si.update(); //distributes/copies the solution of the global on the shared interface > local processor vectors - "repeated vector"
    }

    // Read output folder name from input file
    std::string output_folder_name = data("output_name", "Output_default");
    std::string output_number = data("output_number", "");
    std::string test = output_folder_name + output_number;
    BeatIt::createOutputFolder(test);

    // Output file name
    std::string output_file_name = data("output_file", "/out");


    // case 1 -
    //equation_systems.get_system("Poisson").assemble_before_solve = false;
    // assemble_poisson(parameters);
     // --------------------------




    // Solve the system "Poisson".  Note that calling this
    // member will assemble the linear system and invoke
    // the default numerical solver.  With PETSc the solver can be
    // controlled from the command line.  For example,
    // you can invoke conjugate gradient with:
    //
    // ./introduction_ex3 -ksp_type cg
    //
    // You can also get a nice X-window that monitors the solver
    // convergence with:
    //
    // ./introduction-ex3 -ksp_xmonitor
    //
    // if you linked against the appropriate X libraries when you
    // built PETSc.

    const DofMap& dofmap = system.get_dof_map(); // replicating assembly function

    std::cout << " Looping over elements to set up subdomain IDs"<< std::endl;

    for (const auto& elem: mesh.active_local_element_ptr_range()) { // similar to dealii   - libMesh::Elem* instead of auto&
        int blockid = elem->subdomain_id();//(*elem).subdomain_id();
        libMesh::Point centroid = elem->centroid();
        std::vector<double> u(N);
        std::vector<libMesh::Gradient> du(N);

        // loop over the number of laplace problems we are solving
        for (int i = 0; i < N; i++)
        {
            std::string ui = "u" + std::to_string(i);    //how cool! to use i in the name here!!!
            ExplicitSystem& si = equation_systems.get_system<ExplicitSystem>(ui);
            u[i] = si.point_value(si.variable_number(ui), centroid, *elem);
            du[i] = si.point_gradient(si.variable_number(ui), centroid, *elem);   //may ask for pointers - maybe error
        }

        //std::cout<< "u index= " <<u[tpv1[0]] <<", tpv1[0]= " <<  tpv1[0] <<", tpv1[1]= "<< tpv1[1]<<", tpv1[2]= " << tpv1[2]<<std::endl;
        // setup subregion
        if ((u[tpv1[0]]>tpv1[1])*(u[tpv1[0]]<tpv1[2])){//u[1] > 0.8) { //pv   //(u[tpv1[1]]>tpv1[2])*(u[tpv1[1]]>tpv1[3])
            blockid = 12;
        }
        else if ((u[tpv2[0]]>tpv2[1])*(u[tpv2[0]]<tpv2[2])) { //pv
            blockid = 1;
        }
        else if ((u[tpv3[0]]>tpv3[1])*(u[tpv3[0]]<tpv3[2])) { //pv
            blockid = 2;
        }
        else if ((u[tpv4[0]]>tpv4[1])*(u[tpv4[0]]<tpv4[2])) { //pv
            blockid = 3;
        }
        else if ((u[tfloor[0]]>tfloor[1])*(u[tfloor[0]]<tfloor[2])) { //floor
            blockid = 4;//10;
        }
        else if ((u[tlaa[0]]>tlaa[1])*(u[tlaa[0]]<tlaa[2])) { //laa
            blockid = 5;//20;
        }
        else if ((u[tantra1[0]]>tantra1[1])*(u[tantra1[0]]<tantra1[2])) { //antra
            blockid = 6;//40;
        }
        else if ((u[tantra2[0]]>tantra2[1])*(u[tantra2[0]]<tantra2[2])) { //antra /*
            blockid = 8;//41;
        }
        else if ((u[tseptum[0]]>tseptum[1])*(u[tseptum[0]]<tseptum[2])) { //septum
            blockid = 9;//50;
        }
        else if ((u[tanterior[0]]>tanterior[1])*(u[tanterior[0]]<tanterior[2])) { //anterior
            blockid = 10;//60;
        }
        else if ((u[tposterior[0]]>tposterior[1])*(u[tposterior[0]]<tposterior[2])) { //posterior
            blockid = 11;//70;
        }//*/
        else {               //lateral
            blockid = 7;//12;//80;
        }

//      std::cout <<
        elem->Elem::subdomain_id()=blockid;
        //elem->subdomain_id() = blockid;


        // setup fibers
        libMesh::Gradient f0, s0, n0;
        // SR -> LA
        // Solve the transmural problem as the first one, such that
        // the next line is always true
        s0 = du[0].unit();

        //test
        //std::cout << blockid <<std::endl;

        /*																					UNCOMMENT HERE - THIS IS THE FIRST WAY i DID FIBERS
        switch (blockid)
        {
        //PV left front
        case 12:
        {
				n0 = du[1].unit();
				f0 = s0.cross(n0);
				break;
        }
        //PV left back
        case 1:
        {
				n0 = du[1].unit();
				f0 = s0.cross(n0);
				break;
        }
        //PV right front
        case 2:
        {
				n0 = du[6].unit();
				f0 = s0.cross(n0);
				break;
        }
        //PV right front
        case 3:
        {
				n0 = du[6].unit();
				f0 = s0.cross(n0);
				break;
        }
        // Floor
        case 4:
        {
        	if(u[0]>0.5){
				n0 = du[2].unit();
				f0 = s0.cross(n0);
				break;
        	}
        	else{
				f0 = du[2].unit();
        		n0 = s0.cross(f0);
        		break;
        	}
        }
        // LAA
        case 5:
        {
        	if(u[0]>0.5){
				n0 = du[5].unit();
				f0 = s0.cross(n0);
				break;
        	}
        	else{
				f0 = du[5].unit();
        		n0 = s0.cross(f0);
        		break;
        	}
        }
        //Antra between pv2 and pv3
        case 6:
        {
				n0 = du[6].unit();
				f0 = s0.cross(n0);
				break;
        }
        //Lateral
        case 7:
        {
        	if(u[0]>0.5){
				n0 = du[6].unit();
				f0 = s0.cross(n0);
				break;
        	}
        	else{
				f0 = du[6].unit();
        		n0 = s0.cross(f0);
        		break;
        	}
        }
        //Antra between pv0 and pv1
        case 8:
        {
				n0 = du[1].unit();
				f0 = s0.cross(n0);
				break;
        }
        //Septum
        case 9:
        {
        	if(u[0]>0.5){
				n0 = du[6].unit();
				f0 = s0.cross(n0);
				break;
        	}
        	else{
				f0 = du[6].unit();
        		n0 = s0.cross(f0);
        		break;
        	}
        }
        //Anterior
        case 10:
        {
				n0 = du[6].unit();
				f0 = s0.cross(n0);
				break;
        }
        //Posterior
        case 11:
        {

				n0 = du[6].unit();
				f0 = s0.cross(n0);
				break;
        }
        default:
        {
        	std::cout << "Case default, element = "<< elem <<"\n";
            f0(0) = 1.0; f0(1) = 0.0; f0(2) = 0.0;
            s0(0) = 0.0; s0(1) = 1.0; s0(2) = 0.0;
            n0(0) = 0.0; n0(1) = 0.0; n0(2) = 1.0;
            break;
        }
        }
*/



        switch (blockid)
        {
        //PV left front
        case 12:
        {
				n0 = du[1].unit();
				f0 = s0.cross(n0);
				break;
        }
        //PV left back
        case 1:
        {
				n0 = du[1].unit();
				f0 = s0.cross(n0);
				break;
        }
        //PV right front
        case 2:
        {
				n0 = du[6].unit();
				f0 = s0.cross(n0);
				break;
        }
        //PV right front
        case 3:
        {
				n0 = du[6].unit();
				f0 = s0.cross(n0);
				break;
        }
        // Floor
        case 4:
        {
        	if(u[0]>0.5){
				n0 = du[9].unit();
				f0 = s0.cross(n0);
				break;
        	}
        	else{
				f0 = du[2].unit();
        		n0 = s0.cross(f0);
        		break;
        	}
        }
        // LAA
        case 5:
        {
        	if(u[0]>0.5){
				n0 = du[5].unit();
				f0 = s0.cross(n0);
				break;
        	}
        	else{
				f0 = du[5].unit();
        		n0 = s0.cross(f0);
        		break;
        	}
        }
        //Antra between pv2 and pv3
        case 6:
        {
				n0 = du[6].unit();
				f0 = s0.cross(n0);
				break;
        }
        //Lateral
        case 7:
        {
        	if(u[0]>0.5){
				n0 = du[9].unit();
				f0 = s0.cross(n0);
				break;
        	}
        	else{
				f0 = du[6].unit();
        		n0 = s0.cross(f0);
        		break;
        	}
        }
        //Antra between pv0 and pv1
        case 8:
        {
				n0 = du[1].unit();
				f0 = s0.cross(n0);
				break;
        }
        //Septum
        case 9:
        {
        	if(u[0]>0.5){
				n0 = du[9].unit();
				f0 = s0.cross(n0);
				break;
        	}
        	else{
				f0 = du[6].unit();
        		n0 = s0.cross(f0);
        		break;
        	}
        }
        //Anterior
        case 10:
        {
//				n0 = du[9].unit();
//				f0 = s0.cross(n0);
			f0 = du[8].unit();
			n0 = s0.cross(f0);

				break;
        }
        //Posterior
        case 11:
        {
				n0 = du[8].unit();
				f0 = s0.cross(n0);
				break;
        }
        default:
        {
        	std::cout << "Case default, element = "<< elem <<"\n";
            f0(0) = 1.0; f0(1) = 0.0; f0(2) = 0.0;
            s0(0) = 0.0; s0(1) = 1.0; s0(2) = 0.0;
            n0(0) = 0.0; n0(1) = 0.0; n0(2) = 1.0;
            break;
        }
        }



        /*
        f0.print();
        std::cout<< std::endl;
		*/
        f0 = f0.unit();
        s0=s0.unit();
        n0=n0.unit();

        /*
        std::cout<< "--after .unit\n";
        f0.print();
        std::cout<< std::endl;
		*/

        std::vector<dof_id_type> fibers_dof_indices; //putting this back to the system so that we can export it
        fiber_system.get_dof_map().dof_indices(elem, fibers_dof_indices);
        for (int idim = 0; idim < 3; ++idim)
        {
            fiber_system.solution->set(fibers_dof_indices[idim], f0(idim));
            sheet_system.solution->set(fibers_dof_indices[idim], s0(idim));
            xfiber_system.solution->set(fibers_dof_indices[idim], n0(idim));
        }
    }


    // After solving the system write the solution
    // to a VTK-formatted plot file.
    std::string output_path = test + output_file_name + output_number + ".pvtu";
    VTKIO(mesh).write_equation_systems(output_path, equation_systems);

    // Export fibers in nemesis format
    ExodusII_IO nemesis_exporter(mesh); // creates the exporter
    std::vector<std::string> output_variables(9); // creates a vector that stores the names of vectors
    output_variables[0] = "fibersx";
    output_variables[1] = "fibersy";
    output_variables[2] = "fibersz";
    output_variables[3] = "sheetsx";
    output_variables[4] = "sheetsy";
    output_variables[5] = "sheetsz";
    output_variables[6] = "xfibersx";
    output_variables[7] = "xfibersy";
    output_variables[8] = "xfibersz";
    nemesis_exporter.set_output_variables(output_variables);
    nemesis_exporter.write_equation_systems(test + output_file_name + output_number + ".e", equation_systems); // saves the variables at the nodes
    nemesis_exporter.write_element_data(equation_systems); // this saves the variables at the centroid of the elements
    // All done.
    return 0;
}



// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void assemble_poisson(EquationSystems& es,
    const std::string& libmesh_dbg_var(system_name))
{

    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "Poisson");

    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Poisson");

    //Zero out the solution
    system.solution->zero();

    //Zero out the rhs
    system.rhs->zero();

    //Zero out the solution
    system.matrix->zero();

    // A reference to the  DofMap object for this system.  The  DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the  DofMap
    // in future examples.
    const DofMap& dof_map = system.get_dof_map();

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = dof_map.variable_type(0);                        //??????

    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.  Introduction Example 4
    // describes some advantages of  std::unique_ptr's in the context of
    // quadrature rules.
    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

    // A 5th order Gauss quadrature rule for numerical integration.
    QGauss qrule(dim, FIFTH);

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Declare a special finite element object for
    // boundary integration.
    std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));

    // Boundary integration requires one quadrature rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    QGauss qface(dim - 1, FIFTH);

    // Tell the finite element object to use our
    // quadrature rule.
    fe_face->attach_quadrature_rule(&qface);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe->get_JxW();

    // The physical XY locations of the quadrature points on the element.
    // These might be useful for evaluating spatially varying material
    // properties at the quadrature points.
    const std::vector<Point>& q_point = fe->get_xyz();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>>& phi = fe->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient>>& dphi = fe->get_dphi();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".  These datatypes are templated on
    //  Number, which allows the same code to work for real
    // or complex numbers.
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
    //
    // Element iterators are a nice way to iterate through all the
    // elements, or all the elements that have some property.  The
    // iterator el will iterate from the first to the last element on
    // the local processor.  The iterator end_el tells us when to stop.
    // It is smart to make this one const so that we don't accidentally
    // mess it up!  In case users later modify this program to include
    // refinement, we will be safe and will only consider the active
    // elements; hence we use a variant of the active_elem_iterator.
    for (const auto& elem : mesh.active_local_element_ptr_range())
    {
        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices(elem, dof_indices);

        // Cache the number of degrees of freedom on this element, for
        // use as a loop bound later.  We use cast_int to explicitly
        // convert from size() (which may be 64-bit) to unsigned int
        // (which may be 32-bit but which is definitely enough to count
        // *local* degrees of freedom.
        const unsigned int n_dofs =
            cast_int<unsigned int>(dof_indices.size());

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe->reinit(elem);

        // With one variable, we should have the same number of degrees
        // of freedom as shape functions.
        libmesh_assert_equal_to(n_dofs, phi.size());

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).

        // The  DenseMatrix::resize() and the  DenseVector::resize()
        // members will automatically zero out the matrix  and vector.
        Ke.resize(n_dofs, n_dofs);

        Fe.resize(n_dofs);

        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {

            // Now we will build the element matrix.  This involves
            // a double loop to integrate the test functions (i) against
            // the trial functions (j).
            for (unsigned int i = 0; i != n_dofs; i++)
                for (unsigned int j = 0; j != n_dofs; j++)
                {
                    Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                }

            // This is the end of the matrix summation loop
            // Now we build the element right-hand-side contribution.
            // This involves a single loop in which we integrate the
            // "forcing function" in the PDE against the test functions.
            {
                const Real x = q_point[qp](0);
                const Real y = q_point[qp](1);
                const Real eps = 1.e-3;


                // "fxy" is the forcing function for the Poisson equation.
                // In this case we set fxy to be a finite difference
                // Laplacian approximation to the (known) exact solution.
                //
                // We will use the second-order accurate FD Laplacian
                // approximation, which in 2D is
                //
                // u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
                //                u(i-1,j) + u(i+1,j) +
                //                -4*u(i,j))/h^2
                //
                // Since the value of the forcing function depends only
                // on the location of the quadrature point (q_point[qp])
                // we will compute it here, outside of the i-loop
                const Real fxy = 0;//-(exact_solution(x, y-eps) +
                                 //  exact_solution(x, y+eps) +
                                 //  exact_solution(x-eps, y) +
                                 //  exact_solution(x+eps, y) -
                                 //  4.*exact_solution(x, y))/eps/eps;

                for (unsigned int i = 0; i != n_dofs; i++)
                {
                    Fe(i) += JxW[qp] * fxy * phi[i][qp];
                }
            }
        }

        // We have now reached the end of the RHS summation,
        // and the end of quadrature point loop, so
        // the interior element integration has
        // been completed.  However, we have not yet addressed
        // boundary conditions.  For this example we will only
        // consider simple Dirichlet boundary conditions.
        //
        // There are several ways Dirichlet boundary conditions
        // can be imposed.  A simple approach, which works for
        // interpolary bases like the standard Lagrange polynomials,
        // is to assign function values to the
        // degrees of freedom living on the domain boundary. This
        // works well for interpolary bases, but is more difficult
        // when non-interpolary (e.g Legendre or Hierarchic) bases
        // are used.
        //
        // Dirichlet boundary conditions can also be imposed with a
        // "penalty" method.  In this case essentially the L2 projection
        // of the boundary values are added to the matrix. The
        // projection is multiplied by some large factor so that, in
        // floating point arithmetic, the existing (smaller) entries
        // in the matrix and right-hand-side are effectively ignored.
        //
        // This amounts to adding a term of the form (in latex notation)
        //
        // \frac{1}{\epsilon} \int_{\delta \Omega} \phi_i \phi_j = \frac{1}{\epsilon} \int_{\delta \Omega} u \phi_i
        //
        // where
        //
        // \frac{1}{\epsilon} is the penalty parameter, defined such that \epsilon << 1

        {


            // The following loop is over the sides of the element.
            // If the element has no neighbor on a side then that
            // side MUST live on a boundary of the domain.
            for (auto side : elem->side_index_range())
                if (elem->neighbor_ptr(side) == nullptr)
                {
                    // get boundary info
                    auto boundaryid = mesh.get_boundary_info().boundary_id(elem, side);
                    //std::cout << boundaryid << "\n";

                    //check if boundaryid of element elem is in the Dirichlet BC ID list

                    double elem_dirichlet_bc=0;
                    bool found_side_set_in_list =0;
                    double radius;

                    //Check if side set is on the Idirichlet_side_set_list from input
                    for(int it=0; it<poisson_ptr->dirichlet_side_set_list.size(); it++)
                    {
                    	if(poisson_ptr->dirichlet_side_set_list[it]==boundaryid){
                    	elem_dirichlet_bc=poisson_ptr->dirichlet_bc_list[it];
                    	found_side_set_in_list =1;
                    	radius = poisson_ptr->dirichlet_id_radius[it];
                    	//std::cout<< radius <<std::endl;
                    	}
                    }


                    // The value of the shape functions at the quadrature
                    // points.
                    const std::vector<std::vector<Real>>& phi_face = fe_face->get_phi();

                    // The Jacobian * Quadrature Weight at the quadrature
                    // points on the face.
                    const std::vector<Real>& JxW_face = fe_face->get_JxW();

                    // The XYZ locations (in physical space) of the
                    // quadrature points on the face.  This is where
                    // we will interpolate the boundary value function.
                    const std::vector<Point>& qface_point = fe_face->get_xyz();

                    // Compute the shape function values on the element
                    // face.
                    fe_face->reinit(elem, side);

                    // Some shape functions will be 0 on the face, but for
                    // ease of indexing and generality of code we loop over
                    // them anyway
                    libmesh_assert_equal_to(n_dofs, phi_face.size());

                    // Loop over the face quadrature points for integration.
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        // The location on the boundary of the current
                        // face quadrature point.
                        const Real xf = qface_point[qp](0);
                        const Real yf = qface_point[qp](1);
                        bool is_on_neighborhood0 = 0;
                        bool is_on_neighborhood1 = 0;
                        bool is_on_neighborhood05 = 0;
                        bool is_on_neighborhood = 0;

                        // Verify if the quadrature points are within the neighborhood of one of the id05 list points
                        for (int jj = 0; jj < poisson_ptr->dirichlet_id_points.size(); jj++) {
                            libMesh::Point p_id(qface_point[qp] - poisson_ptr->dirichlet_id_points[jj]);
                            is_on_neighborhood = (p_id.norm() <= poisson_ptr->dirichlet_id_radius[jj]) ? 1 : 0;
                            if (is_on_neighborhood == 1){
                            	elem_dirichlet_bc=poisson_ptr->dirichlet_bc_list[jj];
                                radius = poisson_ptr->dirichlet_id_radius[jj];
                                break;
                            }
                       }

                        // The penalty value.  \frac{1}{\epsilon}
                        // in the discussion above.
                        const Real penalty = 1.e10;

                        // The boundary value.
                        const Real value = exact_solution(xf, yf);


                        // Matrix contribution of the L2 projection.
                        if((found_side_set_in_list && radius<0 ) || (found_side_set_in_list && radius>0 &&  is_on_neighborhood))
                        {
                            for (unsigned int i = 0; i != n_dofs; i++)
                            {
                                for (unsigned int j = 0; j != n_dofs; j++)
                                    Ke(i, j) += JxW_face[qp] * penalty * phi_face[i][qp] * phi_face[j][qp];
                            }
                        }
                        // Right-hand-side contribution of the L2
                        // projection.
                        for (unsigned int i = 0; i != n_dofs; i++) {
                             if((found_side_set_in_list && radius<0 ) || (found_side_set_in_list && radius>0 &&  is_on_neighborhood))
                                Fe(i) += JxW_face[qp] * penalty * elem_dirichlet_bc * phi_face[i][qp]; //dirichlet
                            else
                                Fe(i) += 0;         //neumann
                        }
                    }
                }
        }

        // We have now finished the quadrature point loop,
        // and have therefore applied all the boundary conditions.

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The  SparseMatrix::add_matrix()
        // and  NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    }

    // All done!
}

