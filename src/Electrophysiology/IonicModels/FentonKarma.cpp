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
 * \file FentonKarma.cpp
 *
 * \class FentonKarma
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

#include "Electrophysiology/IonicModels/FentonKarma.hpp"
#include "libmesh/getpot.h"

namespace BeatIt
{


IonicModel* createFentonKarma()
{
    return new FentonKarma;
}

FentonKarma::FentonKarma()
 : super(3, 0, "FentonKarma")
 , M_tau_v_p(3.33)
 , M_tau_v1_m(19.6)
 , M_tau_v2_m(1000.0)
 , M_tau_w_p(667.0)
 , M_tau_w_m(11.0)
 , M_tau_d(0.25)
 , M_tau_0(8.3)
 , M_tau_r(50.0)
 , M_tau_si(45.0)
 , M_kappa(10.0)
 , M_V_c_si(0.85)
 , M_V_c(0.13)
 , M_V_v(0.055)
{
    // Without potential
    M_variablesNames[0] = "v";
    M_variablesNames[1] = "w";
}

void
FentonKarma::setup(GetPot& data, std::string sect)
{
	std::string section = sect + "/FentonKarma";
	// From Fenton: Mechanisms of Spiral break up
    // parameter set 1, 3, 4, 5, 6, 8, 9, 10
	// Parameter set 100 for the Right Atria
	// from: A simple model of the right atrium of the human heart with the sinoatrial and atrioventricular nodes included
	//
	int parameter_set = data(section+"/param_set",       1);
    switch(parameter_set)
    {
        case 1:
        {
            M_tau_v_p = 3.33;
            M_tau_v1_m = 19.6;
            M_tau_v2_m = 1000.0;
            M_tau_w_p = 667.0;
            M_tau_w_m = 11.0;
            M_tau_d = 0.25;
            M_tau_0 = 8.3;
            M_tau_r = 50.0;
            M_tau_si = 45.0;
            M_kappa = 10.0;
            M_V_c_si = 0.85;
            M_V_c = 0.13;
            M_V_v = 0.055;
            break;
        }
        case 3:
        {
            M_tau_v2_m = 1250.0;
            M_tau_w_p = 870.0;
            M_tau_w_m = 41.0;
            M_tau_0 = 12.5;
            M_tau_r = 33.33;
            M_tau_si = 29.0;
            M_V_v = 0.04;
            break;
        }
        case 4:
        {
            M_tau_v_p = 3.33;
            M_tau_v1_m = 15.6;
            M_tau_v2_m = 5.0;
            M_tau_w_p = 350.0;
            M_tau_w_m = 80.0;
            M_tau_d = 0.407;
            M_tau_0 = 9;
            M_tau_r = 34.0;
            M_tau_si = 26.5;
            M_kappa = 15.0;
            M_V_c_si = 0.45;
            M_V_c = 0.15;
            M_V_v = 0.04;
            break;
        }
        case 5:
        {
            M_tau_v_p = 3.33;
            M_tau_v1_m = 12;
            M_tau_v2_m = 2;
            M_tau_w_p = 1000.0;
            M_tau_w_m = 100.0;
            M_tau_d = 0.362;
            M_tau_0 = 5;
            M_tau_r = 33.33;
            M_tau_si = 29.0;
            M_kappa = 15.0;
            M_V_c_si = 0.7;
            M_V_c = 0.13;
            M_V_v = 0.04;
            break;
        }
        case 6:
        {
            M_tau_v_p = 3.33;
            M_tau_v1_m = 9;
            M_tau_v2_m = 8;
            M_tau_w_p = 250.0;
            M_tau_w_m = 60.0;
            M_tau_d = 0.395;
            M_tau_0 = 9;
            M_tau_r = 33.33;
            M_tau_si = 29.0;
            M_kappa = 15.0;
            M_V_c_si = 0.5;
            M_V_c = 0.13;
            M_V_v = 0.04;
            break;
        }
        case 8:
        {
            M_tau_v_p = 13.03;
            M_tau_v1_m = 19.6;
            M_tau_v2_m = 1250;
            M_tau_w_p = 800.0;
            M_tau_w_m = 40.0;
            M_tau_d = 0.45;
            M_tau_0 = 12.5;
            M_tau_r = 33.25;
            M_tau_si = 29.0;
            M_kappa = 10.0;
            M_V_c_si = 0.85;
            M_V_c = 0.13;
            M_V_v = 0.04;
            break;
        }
        case 9:
        {
            M_tau_v_p = 3.33;
            M_tau_v1_m = 15.;
            M_tau_v2_m = 2;
            M_tau_w_p = 670.0;
            M_tau_w_m = 61.0;
            M_tau_d = 0.25;
            M_tau_0 = 12.5;
            M_tau_r = 28.;
            M_tau_si = 29.0;
            M_kappa = 10.0;
            M_V_c_si = 0.45;
            M_V_c = 0.13;
            M_V_v = 0.05;
            break;
        }
        case 10:
        {
            M_tau_v_p = 10.;
            M_tau_v1_m = 40.;
            M_tau_v2_m = 333;
            M_tau_w_p = 1000.0;
            M_tau_w_m = 65.0;
            M_tau_d = 0.115;
            M_tau_0 = 12.5;
            M_tau_r = 25.;
            M_tau_si = 22.22;
            M_kappa = 10.0;
            M_V_c_si = 0.85;
            M_V_c = 0.13;
            M_V_v = 0.025;
            break;
        }
        // Right Atria
        case 100:
        {
            M_tau_v_p = 10.;
            M_tau_v1_m = 18.2;
            M_tau_v2_m = 18.2;
            M_tau_w_p = 1020.0;
            M_tau_w_m = 80.0;
            M_tau_d = 51.0/4.4;
            M_tau_0 = 12.5;
            M_tau_r = 130.;
            M_tau_si = 127;
            M_kappa = 10.0;
            M_V_c_si = 0.85;
            M_V_c = 0.13;
            M_V_v = 0.025;
            break;
        }
        default:
        {
        	std::cout << "Fenton Karma parameter set " << parameter_set << "not found" << std::endl;
        	throw std::runtime_error("Parametr set not found");
        	break;
        }
    }
	double tau_v_p   = data(section+"/tau_v_p",       -1.0);
    double tau_v1_m = data(section+"/tau_v1_m",       -1.0);
    double tau_v2_m = data(section+"/tau_v2_m",       -1.0);
    double tau_w_p = data(section+"/tau_w_p",       -1.0);
    double tau_w_m = data(section+"/tau_w_m",       -1.0);
    double tau_d = data(section+"/tau_d",       -1.0);
    double tau_0 = data(section+"/tau_0",       -1.0);
    double tau_r = data(section+"/tau_r",       -1.0);
    double tau_si = data(section+"/tau_si",       -1.0);
    double kappa = data(section+"/kappa",       -1.0);
    double V_c_si = data(section+"/V_c_si",       -1.0);
    double V_c = data(section+"/V_c",       -1.0);
    double V_v = data(section+"/V_v",       -1.0);
    if(tau_v_p > 0 ) M_tau_v_p = tau_v_p;
    if(tau_v1_m > 0 ) M_tau_v1_m = tau_v1_m;
    if(tau_v2_m > 0 ) M_tau_v2_m = tau_v2_m;
    if(tau_w_p > 0 ) M_tau_w_p = tau_w_p;
    if(tau_w_m > 0 ) M_tau_w_m = tau_w_m;
    if(tau_d > 0 ) M_tau_d = tau_d;
    if(tau_0 > 0 ) M_tau_0 = tau_0;
    if(tau_r > 0 ) M_tau_r = tau_r;
    if(tau_si > 0 ) M_tau_si = tau_si;
    if(kappa > 0 ) M_kappa = kappa;
    if(V_c_si > 0 ) M_V_c_si = V_c_si;
    if(V_c > 0 ) M_V_c = V_c;
    if(V_c > 0 ) M_V_v = V_c;

}


void
FentonKarma::solve(std::vector<double>& variables, double appliedCurrent, double dt)
{
    updateVariables(variables, appliedCurrent, dt);
    variables[0] += dt * evaluateIonicCurrent(variables, appliedCurrent, dt);
}

void
FentonKarma::updateVariables(std::vector<double>& variables, double appliedCurrent, double dt)
{
    double V = variables[0];
    double v = variables[1];
    double w = variables[2];
    double tau_v_m = ( 1 - q(V) ) * M_tau_v1_m + q(V) * M_tau_v2_m;
    variables[1] += dt * ( ( 1-p(V) ) * ( 1 - v ) / tau_v_m - p(V) * v / M_tau_v_p );
    variables[2] += dt * ( ( 1-p(V) ) * ( 1 - w ) / M_tau_w_m - p(V) * w / M_tau_w_p );
}

void
void updateVariables(std::vector<double>& variables, std::vector<double>& rhs, double appliedCurrent, double dt, bool overwrite)
{
    double V = variables[0];
    double v = variables[1];
    double w = variables[2];
    double tau_v_m = ( 1 - q(V) ) * M_tau_v1_m + q(V) * M_tau_v2_m;
    rhs[1] = ( ( 1-p(V) ) * ( 1 - v ) / tau_v_m - p(V) * v / M_tau_v_p );
    rhs[2] = ( ( 1-p(V) ) * ( 1 - w ) / M_tau_w_m - p(V) * w / M_tau_w_p );

  double tau_v_m = ( 1 - q(V) ) * M_tau_v1_m + q(V) * M_tau_v2_m;

  if(overwrite)
  {
      variables[1] += dt * rhs[1];
      variables[2] += dt * rhs[2];

  }

}



void
FentonKarma::updateVariables(double V, std::vector<double>& variables, double dt)
{
	double v = variables[0] ;
    double w = variables[1] ;
    double tau_v_m = ( 1 - q(V) ) * M_tau_v1_m + q(V) * M_tau_v2_m;
    variables[1] += dt * ( ( 1-p(V) ) * ( 1 - v ) / tau_v_m - p(V) * v / M_tau_v_p );
    variables[2] += dt * ( ( 1-p(V) ) * ( 1 - w ) / M_tau_w_m - p(V) * w / M_tau_w_p );
}

double
FentonKarma::evaluateIonicCurrent(std::vector<double>& variables, double appliedCurrent, double dt)
{
    double V = variables[0];
    double v = variables[1];
    double w = variables[2];

    double Ifi = - v * p(V) * (V- M_V_c) * (1-V) / M_tau_d;
    double Iso = V * ( 1 - p(V) ) / M_tau_0 + p(V) / M_tau_r;
    double Isi = - w * ( 1 + std::tanh(M_kappa * ( V - M_V_c_si ) ) ) / 2.0 / M_tau_si;
    return  - ( Ifi + Iso + Isi ) + appliedCurrent;
}
double
FentonKarma::evaluateIonicCurrent(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent, double dt)
{
    double V = v_n[0];
    double v = v_n[1];
    double w = v_n[2];
    double Ifi = - v * p(V) * (V- M_V_c) * (1-V) / M_tau_d;
    double Iso = V * ( 1 - p(V) ) / M_tau_0 + p(V) / M_tau_r;
    double Isi = - w * ( 1 + std::tanh(M_kappa * ( V - M_V_c_si ) ) ) / 2.0 / M_tau_si;

   double f_n =  - ( Ifi + Iso + Isi ) + appliedCurrent;
   V = v_np1[0];
   v = v_np1[1];
   w = v_np1[2];

   Ifi = - v * p(V) * (V- M_V_c) * (1-V) / M_tau_d;
   Iso = V * ( 1 - p(V) ) / M_tau_0 + p(V) / M_tau_r;
   Isi = - w * ( 1 + std::tanh(M_kappa * ( V - M_V_c_si ) ) ) / 2.0 / M_tau_si;
   double f_np1 =  - ( Ifi + Iso + Isi ) + appliedCurrent;

    return 0.5 * (f_n+f_np1);

}

double
FentonKarma::evaluateIonicCurrent(double V, std::vector<double>& variables, double appliedCurrent, double dt)
{
    double v = variables[0];
    double w = variables[1];
    double Ifi = - v * p(V) * (V- M_V_c) * (1-V) / M_tau_d;
    double Iso = V * ( 1 - p(V) ) / M_tau_0 + p(V) / M_tau_r;
    double Isi = - w * ( 1 + std::tanh(M_kappa * ( V - M_V_c_si ) ) ) / 2.0 / M_tau_si;
    return  - ( Ifi + Iso + Isi ) + appliedCurrent;
}

double
FentonKarma::evaluatedIonicCurrent(std::vector<double>& variables, double appliedCurrent, double dt, double h)
{
    double V = variables[0];
    double v = variables[1];
    double w = variables[2];
    double dIfi = - v * p(V) * (1-V) / M_tau_d + v * p(V) * (V- M_V_c) / M_tau_d;
    double dIso = ( 1 - p(V) ) / M_tau_0;
    double tan = std::tanh(M_kappa * ( V - M_V_c_si ) );
    double dIsi = M_kappa*w*(tan * tan - 1.0)/(2*M_tau_si);
    return  - (dIfi + dIso + dIsi);
}

double
FentonKarma::evaluatedIonicCurrent( std::vector<double>& variables,
                                     std::vector<double>& rhs,
                                     double dt,
                                     double h )
{
    double V = variables[0];
    double v = variables[1];
    double w = variables[2];
    double Q = old_variables[0];
    double dv = rhs[1];
    double dw = rhs[2];

    double dIfi = - v * p(V) * (1-V) / M_tau_d + v * p(V) * (V- M_V_c) / M_tau_d;
    double dIso = ( 1 - p(V) ) / M_tau_0;
    double tan = std::tanh(M_kappa * ( V - M_V_c_si ) );
    double dIsi = M_kappa*w*(tan * tan - 1.0)/(2*M_tau_si);

    double dIdV =  - (dIfi + dIso + dIsi);

    double dIfidv = - p(V) * (V- M_V_c) * (1-V) / M_tau_d;
    double dIsidw = - ( 1 + std::tanh(M_kappa * ( V - M_V_c_si ) ) ) / 2.0 / M_tau_si;
    return dIdV * Q + dIfidv * dv + dIsidw * dw;

}

void
FentonKarma::initialize(std::vector<double>& variables)
{
    // Potential V
    variables[0] = 0.0;
    // Recovery Variable v
    variables[1] = 1.0;
    // Recovery Variable w
    variables[2] = 1.0;
}

void
FentonKarma::initializeSaveData(std::ostream& output)
{
    // time -  0
    output << "time V v w";
}


double
FentonKarma::evaluateSAC(double v , double I4f)
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
    return 0.0;
}

} /* namespace BeatIt */
