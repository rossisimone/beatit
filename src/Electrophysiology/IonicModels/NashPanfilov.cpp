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
 * \file NashPanfilov.cpp
 *
 * \class NashPanfilov
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

#include "Electrophysiology/IonicModels/NashPanfilov.hpp"
#include "libmesh/getpot.h"

namespace BeatIt
{


IonicModel* createNashPanfilov()
{
    return new NashPanfilov;
}

NashPanfilov::NashPanfilov()
 : super(2, 0, "NashPanfilov")
 , M_mu1(0.12) // Original value 0.12
 , M_mu2(0.3)
 , M_k(8.0)
 , M_a(0.1) // original value 0.1
 , M_b(0.1)
 , M_epsilon(0.01)
{
    // Without potential
    M_variablesNames[0] = "r";
}

void
NashPanfilov::setup(GetPot& data, std::string sect)
{
	std::string section = sect + "/NashPanfilov";
	M_mu1     = data(section+"/mu1",       0.12);
	M_mu2     = data(section+"/mu2",       0.3);
	M_k       = data(section+"/k",         8.0);
	M_a       = data(section+"/a",         0.1);
	M_b       = data(section+"/b",         0.1);
	M_epsilon = data(section+"/epsilon",   0.01);
}


void
NashPanfilov::solve(std::vector<double>& variables, double appliedCurrent, double dt)
{
    updateVariables(variables, appliedCurrent, dt);
    variables[0] += dt * evaluateIonicCurrent(variables, appliedCurrent, dt);
}

void
NashPanfilov::updateVariables(std::vector<double>& variables, double appliedCurrent, double dt)
{
    double V = variables[0];
    double r = variables[1];
    variables[1] += dt * (M_epsilon + M_mu1 * r / (M_mu2 + V) ) *
                         (- r - M_k * V * (V - M_b - 1.0) );
}

void
NashPanfilov::updateVariables(std::vector<double>& variables, std::vector<double>& rhs, double appliedCurrent, double dt, bool overwrite)
{
    double V = v_n[0];
    double r = v_n[1];
    double f_n =  (M_epsilon + M_mu1 * r / (M_mu2 + V) ) *
                             (- r - M_k * V * (V - M_b - 1.0) );
    // we store in rhs[0] Q^n
    rhs[1] = f_n;
    if(overwrite)
    {
        variables[1] += dt*f_n;
    }
}



void
NashPanfilov::updateVariables(double V, std::vector<double>& variables, double dt)
{
	double r = variables[0] ;
    variables[0] += dt * (M_epsilon + M_mu1 * r / (M_mu2 + V) ) *
                         (- r - M_k * V * (V - M_b - 1.0) );
}

double
NashPanfilov::evaluateIonicCurrent(std::vector<double>& variables, double appliedCurrent, double dt)
{
    double V = variables[0];
    double r = variables[1];
    return  - M_k * V * (V - 1.0) * (V - M_a) - r * V + appliedCurrent;
}
double
NashPanfilov::evaluateIonicCurrent(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent, double dt)
{
    double V = v_n[0];
    double r = v_n[1];
   double f_n =  - M_k * V * (V - 1.0) * (V - M_a) - r * V + appliedCurrent;
   V = v_np1[0];
   r = v_np1[1];

   double f_np1 = - M_k * V * (V - 1.0) * (V - M_a) - r * V + appliedCurrent;

    return 0.5 * (f_n+f_np1);

}

double
NashPanfilov::evaluateIonicCurrent(double V, std::vector<double>& variables, double appliedCurrent, double dt)
{
    double r = variables[0];
    return  - M_k * V * (V - 1.0) * (V - M_a) - r * V + appliedCurrent;
}

double
NashPanfilov::evaluatedIonicCurrent(std::vector<double>& variables, double appliedCurrent, double dt, double h)
{
    double V = variables[0];
    double r = variables[1];
    return  - M_k * ( (V - 1.0) * (V - M_a) +  V * (V - M_a) + V * (V - 1.0) ) - r;
}
double
NashPanfilov::evaluatedIonicCurrent( std::vector<double>& variables,
                                     std::vector<double>& gating_rhs,
                                     double dt,
                                     double h )
{
    double V = variables[0];
    double r = variables[1];
    double dIdV =  - M_k * ( (V - 1.0) * (V - M_a) +  V * (V - M_a) + V * (V - 1.0) ) - r;
    double dIdr =  - V;
    double Q = gating_rhs[0];
    double dr = gating_rhs[1];
    return dIdV * Q + dIdr * dr;
}



void
NashPanfilov::initialize(std::vector<double>& variables)
{
    // Potential V
    variables[0] = 0.0;
    // Recovery Variable
    variables[1] = 0.0;
}

void
NashPanfilov::initializeSaveData(std::ostream& output)
{
    // time -  0
    output << "time v r";
}


double
NashPanfilov::evaluateSAC(double v , double I4f)
{
    // From
    // An electromechanical model of cardiac tissue: Constitutive issues and electrophysiological effects
    // C. Cherubinia, S. Filippia, P. Nardinocchib, L. Teresi
    double sac = 0.0;
    double Gs = 5;
    double delta = 60;
    double l_ref = 1.1;
    double l = std::sqrt(I4f);
    double Vt = 0.56;
    sac = Gs * (v - Vt) / ( 1  + std::exp(-delta * (l-l_ref) ) );
    return sac;
}

} /* namespace BeatIt */
