/*
 * Noise.cpp
 *
 *  Created on: Mar 16, 2021
 *      Author: srossi
 */

#include "Noise.hpp"
#include "libmesh/mesh_generation.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/numeric_vector.h"


Noise::Noise(  )
: Nf(1)
, gamma_noise(0.5)
, s(15)
, L_fibers(10)
, PHI_fibers(10*M_PI)
, sigma(1.0)
, gamma(1.0)
, beta(0.0)
, kappa(1.0)
, L(4.0)
, Ly(4.0)
, g_system_ptr (nullptr)
, f_system_ptr (nullptr)
, locator_ptr (nullptr)
, grid_min_x(0)
, grid_min_y(0)
, grid_min_z(0)
, grid_max_x(1)
, grid_max_y(1)
, grid_max_z(1)
, grid_el_x(4)
, grid_el_y(4)
, grid_el_z(0)
, use_white_noise(false)
, use_fiber_noise(false)
, rd()
, mt(rd())
, uniform_dist(0.0, 2*M_PI)
, gauss_dist(0.0, 1.0)
, lattice(nullptr)
, noiseType(NoiseType::PERLIN)
, projection(false)
{
    // TODO Auto-generated constructor stub
}

Noise::~Noise()
{
    if(locator_ptr) delete locator_ptr;
    if(lattice) delete lattice;
    // TODO Auto-generated destructor stub
}

void
Noise::create_lattice(libMesh::EquationSystems& es)
{
    std::cout << "* Noise: create_lattice " << std::endl;

    lattice = new libMesh::ReplicatedMesh(es.comm());
    if(grid_el_z <= 0)
    {
        libMesh::MeshTools::Generation::build_square (*lattice,
                grid_el_x, grid_el_y,
                grid_min_x, grid_max_x,
                grid_min_y, grid_max_y,
                libMesh::QUAD4);
    }
    else
    {
        libMesh::MeshTools::Generation::build_cube (*lattice,
                grid_el_x, grid_el_y, grid_el_z ,
                grid_min_x, grid_max_x,
                grid_min_y, grid_max_y,
                grid_min_z, grid_max_z,
                libMesh::HEX8);
    }
    lattice->print_info();
    std::cout << "* Noise: PointLocatorTree  " << std::endl;
    locator_ptr = new libMesh::PointLocatorTree(*lattice);

}

void
Noise::create_lattice( libMesh::EquationSystems& es,
                        int elX, int elY, int elZ,
                       double xmin, double xmax,
                       double ymin, double ymax,
                       double zmin, double zmax)
{
    grid_el_x = elX;
    grid_el_y = elY;
    grid_el_z = elZ;
    grid_min_x= xmin;
    grid_max_x= xmax;
    grid_min_y=ymin;
    grid_max_y=ymax;
    grid_min_z=zmin;
    grid_max_z=zmax;
    create_lattice(es);
}

void
Noise::setup(GetPot& data, std::string section)
{
    Nf = data(section+"Nf",1);
    gamma_noise = data(section+"gamma_noise", 0.5);
    s = data(section+"s", 15);
    L_fibers = data(section+"L_fibers", 10);
    PHI_fibers = data(section+"PHI_fibers", 10*M_PI);
    sigma = data(section+"sigma", 1.0);
    gamma = data(section+"gamma", 1.0);
    beta = data(section+"beta",0.0);
    kappa = data(section+"kappa",1.0);
    L = data(section+"L",4.0);
    Ly = data(section+"Ly",4.0);
    grid_min_x = data(section+"grid_min_x",0);
    grid_min_y = data(section+"grid_min_y",0);
    grid_min_z = data(section+"grid_min_z",0);
    grid_max_x = data(section+"grid_max_x",1);
    grid_max_y = data(section+"grid_max_y",1);
    grid_max_z = data(section+"grid_max_z",1);
    grid_el_x = data(section+"grid_el_x",4);
    grid_el_y = data(section+"grid_el_y",4);
    grid_el_z = data(section+"grid_el_z",0);
    use_white_noise = data(section+"use_white_noise",false);
    use_fiber_noise = data(section+"use_fiber_noise",false);
}


void
Noise::generate_noise_field(libMesh::EquationSystems& es)
{
    std::cout << "* Noise: generate_noise_field " << std::endl;
    create_lattice(es);
    std::cout << "* Noise: create gradients " << std::endl;
    libMesh::EquationSystems lattice_equation_systems(*lattice);
    libMesh::ExplicitSystem & g_system = lattice_equation_systems.add_system<libMesh::ExplicitSystem>("Perlin");
    g_system_ptr = &g_system;
    // Add the variables "u" & "v" to "Stokes".  They
    // will be approximated using FIRST-order approximation.
    g_system.add_variable("ux", libMesh::FIRST);
    g_system.add_variable("uy", libMesh::FIRST);
    g_system.add_variable("uz", libMesh::FIRST);
    g_system.add_variable("value", libMesh::FIRST);
    lattice_equation_systems.init();
    std::vector<libMesh::dof_id_type> dof_indices;
    // define lattice gradients;
    for (auto & node : lattice->local_node_ptr_range())
    {
        g_system.get_dof_map().dof_indices(node, dof_indices);
        double angle1 = uniform_dist(mt);
        double angle2 = 0.0;
        if(grid_el_z > 0 ) angle2 = uniform_dist(mt);

        double ux = std::cos(angle1) * std::cos(angle2);
        double uy = std::sin(angle1) * std::cos(angle2);
        double uz = std::sin(angle2);

        libMesh::Point u(ux,uy, uz);
        double unorm = u.norm();
        u = u / u.norm();

        g_system.solution->set(dof_indices[0], u(0));
        g_system.solution->set(dof_indices[1], u(1));
        g_system.solution->set(dof_indices[2], u(2));
        g_system.solution->set(dof_indices[3], gauss_dist(mt));
    }

    // Declare the Poisson system and its variables.
    // The Poisson system is another example of a steady system.
    es.add_system<libMesh::ExplicitSystem> ("Noise");
    es.get_system("Noise").add_variable("noise", libMesh::CONSTANT, libMesh::MONOMIAL);
    es.get_system("Noise").init();

    if(projection)
    {
        //es.add_system<LinearImplicitSystem> ("Poisson");
        // Adds the variable "u" to "Poisson".  "u"
        // will be approximated using second-order approximation.
        //es.get_system("Poisson").add_variable("u", FIRST);
        //es.get_system("Poisson").init();
        //assemble_noise_system(es);
        //es.get_system("Poisson").solve();
    }
    else
    {
        evaluateNoise(es);
    }
    es.get_system("Noise").solution->close();
    double u_max =   es.get_system("Noise").solution->max();
    double u_min =   es.get_system("Noise").solution->min();
    double du = u_max - u_min;
    es.get_system("Noise").solution->add(-u_min);
    (*es.get_system("Noise").solution) /= du;


}


void
Noise::evaluateNoise(libMesh::EquationSystems& es)
{
    std::cout << "* Noise: evaluateNoise " << std::endl;

    // find size
    double hx = (grid_max_x - grid_min_x) / grid_el_x;
    double hy = (grid_max_y - grid_min_y) / grid_el_y;
    double hz = (grid_max_z - grid_min_z) / grid_el_z;


    libMesh::ExplicitSystem& system = es.get_system<libMesh::ExplicitSystem> ("Noise");
    const libMesh::DofMap & dof_map = system.get_dof_map();
    std::vector<libMesh::dof_id_type> dof_indices;
    std::vector<libMesh::dof_id_type> g_dof_indices;

    const libMesh::MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    std::vector<double> g_local;
    g_system_ptr->solution->localize(g_local);

    for (const auto & elem : mesh.active_local_element_ptr_range())
    {
        dof_map.dof_indices (elem, dof_indices);

        double noise = 0.0;
        // white noise
        switch(noiseType)
        {
            case NoiseType::PERLIN:
            {
                libMesh::Point x = elem->centroid();
                const libMesh::Elem* mapped_element = (*locator_ptr)(x);
                noise = 0.5;
                if(mapped_element)
                {
                    for( auto & node : mapped_element->node_ref_range())
                    {
                        g_system_ptr->get_dof_map().dof_indices(&node, g_dof_indices);
                        libMesh::Point xi(node);

                       // g_system_ptr->solution->localize(g_local,g_dof_indices);
//                        libMesh::Point gi( g_local[0], g_local[1], g_local[2] );
//                        libMesh::Point gi( (*g_system_ptr->solution)(g_dof_indices[0]),
//                                          (*g_system_ptr->solution)(g_dof_indices[1]),
//                                          (*g_system_ptr->solution)(g_dof_indices[2]) );
                        libMesh::Point gi( g_local[g_dof_indices[0]], g_local[g_dof_indices[1]], g_local[g_dof_indices[2]]);
                        double gc = 1.0;
                        double PNf = 0.0;
                        for(int cf = 0; cf < Nf; ++cf)
                        {
                            libMesh::Point ri(x-xi);
                            double dX = 1-std::abs(ri(0)) / hx;
                            double dY = 1-std::abs(ri(1)) / hy;
                            double dZ = 1-std::abs(ri(2)) / hz;

                            double aX = ((6*dX - 15)*dX + 10)*dX*dX*dX;
                            double aY = ((6*dY - 15)*dY + 10)*dY*dY*dY;
                            double aZ = ((6*dZ - 15)*dZ + 10)*dZ*dZ*dZ;


                            //if(grid_el_z <= 0) P += std::sqrt(2.0) * 0.5 * aX * aY * gi.contract(ri);
                            //else P += 1.0 / std::sqrt(3.0) * aX * aY * aZ* gi.contract(ri);
                            if(grid_el_z <= 0) PNf += gc * aX * aY * gi.contract(ri);
                            else PNf += gc * aX * aY * aZ* gi.contract(ri);
                            gc *= gamma_noise;
                        }
                        noise += (1.0 - gamma_noise) / (1.0 - std::pow(gamma_noise, Nf) ) * ( 0.5 * std::sqrt(dim)) * PNf;
                    }
                }
                break;
            }
            case NoiseType::WHITE:
            default:
            {
                noise =  gauss_dist(mt);
                break;
            }
        }
        system.solution->set(dof_indices[0], noise);

    }
    std::cout << "* Noise: evaluateNoise done " << std::endl;


}


void
Noise::assemble_noise_system(libMesh::EquationSystems& es)
{
//    // Get a constant reference to the mesh object.
//    const MeshBase & mesh = es.get_mesh();
//    // The dimension that we are running
//    const unsigned int dim = mesh.mesh_dimension();
//
//    // Get a reference to the LinearImplicitSystem we are solving
//    LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("Poisson");
//    LinearImplicitSystem & w_system = es.get_system<LinearImplicitSystem> ("Noise");
//
//    // A reference to the  DofMap object for this system.  The  DofMap
//    // object handles the index translation from node and element numbers
//    // to degree of freedom numbers.  We will talk more about the  DofMap
//    // in future examples.
//    const DofMap & dof_map = system.get_dof_map();
//    const DofMap & w_dof_map = w_system.get_dof_map();
//    // Get a constant reference to the Finite Element type
//    // for the first (and only) variable in the system.
//    FEType fe_type = dof_map.variable_type(0);


}
