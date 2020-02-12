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
 * \file Cubic.cpp
 *
 * \class Cubic
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
 * Created on: Aug 11, 2016
 *
 */

#include "Electrophysiology/IonicModels/Cubic.hpp"
#include "libmesh/getpot.h"
#include "math.h"
namespace BeatIt
{


IonicModel* createCubic()
{
    return new Cubic;
}

Cubic::Cubic()
 : super(1, 0, "Cubic")
 , M_alpha(0.5)
 , M_v0(0.0)
, M_v1(0.1)
, M_v2(1.0)
, M_k(1)
{
}

void
Cubic::setup(GetPot& data, std::string sect)
{
    super::setup(data, sect);
	std::string section = sect + "/cubic";
	M_dimensional       = data(section+"/dimensional", false);
	M_v0       = data(section+"/v0",       -85.0);
    M_v1       = data(section+"/v1",       -57.6);
    M_v2       = data(section+"/v2",       30.);
    M_k       = data(section+"/k",         1.4e-3);
	std::cout << "Cubic: " << std::endl;
	std::cout << "alpha: " << M_alpha << std::endl;
	std::cout << "v0: " << M_v0  << std::endl;
    std::cout << "v1: " << M_v1  << std::endl;
    std::cout << "v2: " << M_v2  << std::endl;
    std::cout << " k: " << M_k  << std::endl;
}

void
Cubic::updateVariables(std::vector<double>& variables, double appliedCurrent, double dt)
{
}

void
Cubic::updateVariables(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent, double dt)
{
}

void
Cubic::updateVariables(std::vector<double>& v_n, std::vector<double>& rhs, double appliedCurrent, double dt, bool overwrite)
{
}

void
Cubic::updateVariables(double V, std::vector<double>& variables, double dt)
{
}

double
Cubic::evaluateIonicCurrent(std::vector<double>& variables, double appliedCurrent, double dt)
{
    double V = variables[0];
    return  M_k * (V-M_v0) * (V-M_v1) * (V-M_v2);
}

double
Cubic::evaluateIonicCurrentH(std::vector<double>& variables, double appliedCurrent, double dt, double h)
{
    return  0;
}


double
Cubic::evaluateIonicCurrentTimeDerivative( std::vector<double>& variables,
                                     std::vector<double>& rhs,
                                     double dt,
                                     double h )
{
    double V = variables[0];
    double Q = rhs[0];
    double dIdV =  M_k * ( (V-M_v0) * (V-M_v1) +
                           (V-M_v1) * (V-M_v2) +
                           (V-M_v2) * (V-M_v0) ) ;
    return  dIdV*Q;

}


double
Cubic::evaluateIonicCurrent(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent, double dt)
{
    double V = v_n[0];
       double f_n = M_k * (V-M_v0) * (V-M_v1) * (V-M_v2);
   V = v_np1[0];
      double f_np1 = M_k * (V-M_v0) * (V-M_v1) * (V-M_v2);
    return 0.5 * (f_n+f_np1);

}

double
Cubic::evaluateIonicCurrent(double V, std::vector<double>& variables, double appliedCurrent, double dt)
{
    return   M_k * (V-M_v0) * (V-M_v1) * (V-M_v2);
}


void
Cubic::initialize(std::vector<double>& variables)
{
    // Potential V
    variables[0] = M_v0;
}

void
Cubic::initializeSaveData(std::ostream& output)
{
    // time -  0
    output << "time v";
}



} /* namespace BeatIt */
