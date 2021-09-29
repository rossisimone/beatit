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

#include "Electrophysiology/IonicModels/Lombardo2016.hpp"
#include "libmesh/getpot.h"

namespace BeatIt
{


IonicModel* createLombardo2016()
{
    return new Lombardo2016;
}


Lombardo2016::Lombardo2016()
 : super(4, 0, "Lombardo2016")
,  M_uc(0.13)
,  M_uv(0.04)
,  M_uw(0.13)
,  M_ud(0.13)
,  M_tvm(19.6) // ms
,  M_tvp(3.33) // ms
,  M_twm(41) // ms
,  M_twp(870) // ms
,  M_tsp(1) // ms
,  M_tsm(1) // ms
,  M_ucsi(0.85)
,  M_xk(10)
,  M_td(0.25)  // ms
,  M_to(12.5)  // ms
,  M_tsoa(33.3) // ms
,  M_tsob(33.3) // ms
,  M_uso(0.85)
,  M_xtso(10)  // ms
,  M_tsi(29)
,  M_tvmm(1250)
,  D(1) // ms
,  M_u0(0)
,  M_um(1)
,  M_una(0.23)
{
    // Without potential
    M_variablesNames[0] = "v";
    M_variablesNames[1] = "w";
    M_variablesNames[2] = "d";
}

void
Lombardo2016::setup(GetPot& data, std::string sect)
{
	std::string section = sect + "/Lombardo2016";
	super::setup(data,sect);

    double uc = data(section+"/uc", -1.0);
    double uv = data(section+"/uv", -1.0);
    double uw = data(section+"/uw", -1.0);
    double ud = data(section+"/ud", -1.0);
    double tvm = data(section+"/tvm", -1.0); // ms
	double tvp = data(section+"/tvp", -1.0); // ms
    double twm = data(section+"/twm", -1.0); // ms
	double twp = data(section+"/twp", -1.0); // ms
    double tsm = data(section+"/tsm", -1.0); // ms
    double tsp = data(section+"/tsp", -1.0); // ms
    double ucsi = data(section+"/ucsi", -1.0);
    double xk = data(section+"/xk", -1.0);
	double td = data(section+"/td", -1.0);  // ms
	double to = data(section+"/to", -1.0);  // ms
    double tsoa = data(section+"/tsoa", -1.0);// ms
    double tsob = data(section+"/tsob", -1.0);// ms
    double uso = data(section+"/uso", -1.0);
    double xtso = data(section+"/xtso", -1.0); // ms
    double tsi = data(section+"/tsi", -1.0); // ms
    double tvmm = data(section+"/tvmm", -1.0); // ms
    double DD = data(section+"/D", -1.0); // ms


    M_uc  = ( uc  > 0) ? uc  : M_uc;
    M_uv  = ( uv  > 0) ? uv  : M_uv;
    M_uw  = ( uw  > 0) ? uw  : M_uc;
    M_ud  = ( ud  > 0) ? ud  : M_ud;

    M_tvm  = ( tvm  > 0) ? tvm  : M_tvm;
    M_tvp  = ( tvp  > 0) ? tvp  : M_tvp;
    M_twm  = ( twm  > 0) ? twm  : M_twm;
    M_twp  = ( twp  > 0) ? twp  : M_twp;
    M_tsm  = ( tsm  > 0) ? tsm  : M_tsm;
    M_tsp  = ( tsp  > 0) ? tsp  : M_tsp;

    M_td  = ( td  > 0) ? td  : M_td;
    M_to  = ( to  > 0) ? to  : M_to;

    M_ucsi= ( ucsi> 0) ? ucsi: M_ucsi;

    M_xk= ( xk> 0) ? xk: M_xk;
    M_td = ( td > 0) ? td : M_td;
    M_to = ( to > 0) ? to : M_to;
    M_tsoa = ( tsoa > 0) ? tsoa : M_tsoa;
    M_tsob = ( tsob > 0) ? tsob : M_tsob;

    M_uso = ( uso > 0) ? uso : M_uso;
    M_xtso = ( xtso > 0) ? xtso : M_xtso;
    M_tsi = ( tsi > 0) ? tsi : M_tsi;
    M_tvmm = ( tvmm > 0) ? tvmm : M_tvmm;
    D = ( DD > 0) ? DD : D;

}


void
Lombardo2016::updateVariables(std::vector<double>& variables, double appliedCurrent, double dt)
{
    double u = variables[0];
    double v = variables[1];
    double w = variables[2];
    double d = variables[3];

    double Humuna = H(u-M_una);
    double Humuv = H(u-M_uv);
    double Humuw = H(u-M_uw);
    double Humud = H(u-M_ud);

    double dv = (1 - Humuna) * ( 1 - v ) / ( ( 1 - Humuv) * M_tvm + Humuv * M_tvmm  ) - (Humuna) * v / M_tvp;
    double dw = (1 - Humuw) * ( 1 - w ) / M_twm - (Humuw) * w / M_twp;
    double dd = ( (1 - Humud) / M_tsm + Humud / M_tsp ) *  ( 0.5 * ( 1 + std::tanh( (u-M_ucsi) * M_xk ) ) - d );

    // update v
    variables[1] += dt * dv;
    // update w
    variables[2] += dt * dw;
    // update d
    variables[3] += dt * dd;
}

void
Lombardo2016::updateVariables(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent, double dt)
{
   std::cout << "* Lombardo2016::updateVariables -- NOT IMPLEMENTED" << std::endl;
   throw std::runtime_error("not implemented");
}


void 
Lombardo2016::updateVariables(std::vector<double>& variables, std::vector<double>& rhs, double appliedCurrent, double dt, bool overwrite)
{
    std::cout << "* Lombardo2016::updateVariables -- NOT IMPLEMENTED" << std::endl;
    throw std::runtime_error("not implemented");
}



void
Lombardo2016::updateVariables(double V, std::vector<double>& variables, double dt)
{
    std::cout << "* Lombardo2016::updateVariables -- NOT IMPLEMENTED" << std::endl;
    throw std::runtime_error("not implemented");
}

double
Lombardo2016::evaluateIonicCurrent(std::vector<double>& variables, double appliedCurrent, double dt)
{
    double u = variables[0];
    double v = variables[1];
    double w = variables[2];
    double d = variables[3];

    double Humuc = H(u-M_uc);
    double Humuso = H(u-M_uso);
    double Humuna = H(u-M_una);

    double Ifi = - v * Humuna * ( u - M_una) * (M_um - u ) / M_td;

    double tso = M_tsoa + 0.5 * (M_tsob - M_tsoa) * ( 1 + std::tanh( (u - M_uso) * M_xtso ) );

    double Iso = (u - M_u0) * (1 - Humuc) / M_to + Humuc / tso;

    double Isi = - w * d / M_tsi;
    return  Ifi + Iso + Isi;
}


double
Lombardo2016::evaluateIonicCurrent(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent, double dt)
{
    std::cout << "* Lombardo2016::evaluateIonicCurrent -- NOT IMPLEMENTED" << std::endl;
    throw std::runtime_error("not implemented");

}

double
Lombardo2016::evaluateIonicCurrent(double V, std::vector<double>& variables, double appliedCurrent, double dt)
{
    std::cout << "* Lombardo2016::evaluateIonicCurrent -- NOT IMPLEMENTED" << std::endl;
    throw std::runtime_error("not implemented");
}

double
Lombardo2016::evaluateIonicCurrentTimeDerivative( std::vector<double>& variables,
                                     std::vector<double>& rhs,
                                     double dt,
                                     double h )
{
   return 0;
}

void
Lombardo2016::initialize(std::vector<double>& variables)
{
    // Potential V
    variables[0] = 0.0;
    // Recovery Variable v
    variables[1] = 1.0;
    // Recovery Variable w
    variables[2] = 1.0;
    // Recovery Variable s
    variables[3] = 0.0;
}

void
Lombardo2016::initializeSaveData(std::ostream& output)
{
    // time -  0
    output << "time V v w d\n";
}


double
Lombardo2016::evaluateSAC(double v , double I4f)
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
