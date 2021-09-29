/*
 * Noise.hpp
 *
 *  Created on: Mar 16, 2021
 *      Author: srossi
 */

#ifndef SRC_UTIL_NOISE_HPP_
#define SRC_UTIL_NOISE_HPP_

#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>
#include "libmesh/getpot.h"
#include "libmesh/explicit_system.h"
#include <libmesh/point_locator_tree.h>
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/replicated_mesh.h"

enum class NoiseType { WHITE, PERLIN, VALUE };

class Noise
{
public:
    Noise();
    virtual ~Noise();

    void setup(GetPot& data, std::string section = "");

    void generate_noise_field(libMesh::EquationSystems& es);
protected:

    void evaluateNoise(libMesh::EquationSystems& es);
    void assemble_noise_system(libMesh::EquationSystems& es);
    void create_lattice(libMesh::EquationSystems& es);
    void create_lattice( libMesh::EquationSystems& es,
                         int elX, int elY, int elZ,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         double zmin, double zmax);


    int Nf;
    double gamma_noise;
    double s;
    double L_fibers;
    double PHI_fibers;

    double sigma;
    double gamma;
    double beta;
    double kappa;
    double L;
    double Ly;
    libMesh::ExplicitSystem * g_system_ptr;
    libMesh::ExplicitSystem * f_system_ptr;
    libMesh::PointLocatorTree *  locator_ptr;

    double grid_min_x;
    double grid_min_y;
    double grid_min_z;
    double grid_max_x;
    double grid_max_y;
    double grid_max_z;
    int grid_el_x;
    int grid_el_y;
    int grid_el_z;

    bool use_white_noise;
    bool use_fiber_noise;

    std::random_device rd;
    std::mt19937 mt;
    std::uniform_real_distribution<double> uniform_dist;
    std::normal_distribution<double> gauss_dist;

    libMesh::ReplicatedMesh * lattice;

    NoiseType noiseType;
    bool projection;

};

#endif /* SRC_UTIL_NOISE_HPP_ */
