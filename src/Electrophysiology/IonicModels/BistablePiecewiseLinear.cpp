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
 * \file BistablePiecewiseLinear.cpp
 *
 * \class BistablePiecewiseLinear
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

#include "Electrophysiology/IonicModels/BistablePiecewiseLinear.hpp"
#include "libmesh/getpot.h"
#include "math.h"
namespace BeatIt
{


IonicModel* createBistablePiecewiseLinear()
{
    return new BistablePiecewiseLinear;
}

BistablePiecewiseLinear::BistablePiecewiseLinear()
 : super(1, 0, "BistablePiecewiseLinear")
 , M_alpha(0.5)
 , M_v0(0.5)
{
}

void
BistablePiecewiseLinear::setup(GetPot& data, std::string sect)
{
    super::setup(data, sect);
	std::string section = sect + "/BistablePiecewiseLinear";
	M_alpha       = data(section+"/a",         0.5);
	M_v0       = data(section+"/v0",         0.0);
	std::cout << "BistablePiecewiseLinear: " << std::endl;
	std::cout << "alpha: " << M_alpha << std::endl;
	std::cout << "v0: " << M_v0  << std::endl;
}

void
BistablePiecewiseLinear::updateVariables(std::vector<double>& variables, double appliedCurrent, double dt)
{
}

void
BistablePiecewiseLinear::updateVariables(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent, double dt)
{
}

void
BistablePiecewiseLinear::updateVariables(std::vector<double>& v_n, std::vector<double>& rhs, double appliedCurrent, double dt, bool overwrite)
{
}

void
BistablePiecewiseLinear::updateVariables(double V, std::vector<double>& variables, double dt)
{
}

double
BistablePiecewiseLinear::evaluateIonicCurrent(std::vector<double>& variables, double appliedCurrent, double dt)
{
    double V = variables[0];
    return  V - BistablePiecewiseLinear::Heaviside(V);
//    return  - V + (0.5 + 0.5 * std::tanh(100*(V - M_alpha) ) ) + appliedCurrent;
}

double
BistablePiecewiseLinear::evaluateIonicCurrentH(std::vector<double>& variables, double appliedCurrent, double dt, double h)
{
    double V = variables[0];
    double H = 0.0;
    double xi = V-M_alpha;
    if(xi >= h ) H = 1.0;
    else if( xi > -h && xi < h)
    {
    	H = 0.5 * (xi / h + 1);
    }
    return  V - H;
}


double
BistablePiecewiseLinear::evaluateIonicCurrentTimeDerivative( std::vector<double>& variables,
                                     std::vector<double>& rhs,
                                     double dt,
                                     double h )
{
    double Q = rhs[0];
    double dIdV =  1.0;

    return  dIdV*Q;

}


double
BistablePiecewiseLinear::evaluateIonicCurrent(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent, double dt)
{
    double V = v_n[0];
//   double f_n =- V + BistablePiecewiseLinear::Heaviside(V) + appliedCurrent;
       double f_n = V -  (0.5 + 0.5 * std::tanh(100*(V - M_alpha) ) );
   V = v_np1[0];
//   double f_np1 =- V + BistablePiecewiseLinear::Heaviside(V) + appliedCurrent;
      double f_np1 = V -  (0.5 + 0.5 * std::tanh(100*(V - M_alpha) ) );
    return 0.5 * (f_n+f_np1);

}

double
BistablePiecewiseLinear::evaluateIonicCurrent(double V, std::vector<double>& variables, double appliedCurrent, double dt)
{
    return  V - BistablePiecewiseLinear::Heaviside(V);
}


void
BistablePiecewiseLinear::initialize(std::vector<double>& variables)
{
    // Potential V
    variables[0] = M_v0;
}

void
BistablePiecewiseLinear::initializeSaveData(std::ostream& output)
{
    // time -  0
    output << "time v";
}



} /* namespace BeatIt */
