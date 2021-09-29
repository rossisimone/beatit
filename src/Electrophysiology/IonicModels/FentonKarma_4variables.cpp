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

#include "Electrophysiology/IonicModels/FentonKarma_4variables.hpp"
#include "libmesh/getpot.h"

namespace BeatIt
{


IonicModel* createFentonKarma4v()
{
    return new FentonKarma4v;
}

FentonKarma4v::FentonKarma4v()
 : super(4, 0, "FentonKarma4v")
, M_tvp(3.33) // ms
, M_tv1m(19.2) // ms
, M_tv2m(10.0) // ms
, M_twp(160.0) // ms
, M_tw1m(75.0) // ms
, M_tw2m(75.0) // ms
, M_td(0.065)  // ms
, M_tsi(31.8364) // ms
, M_to(39.0)  // ms
, M_ta(0.009)  // ms
, M_uc(0.23)
, M_uv(0.055)
, M_uw(0.146)
, M_uo(0)
, M_um(1.0)
, M_ucsi(0.8)
, M_uso(0.3)
, M_rsp(0.02)
, M_rsm(1.2)
, M_k(3.0)
, M_aso(0.115)// ms
, M_bso(0.84)// ms
, M_cso(0.02)// ms
, M_pulmonary_vein_model(false)
{
    // Without potential
    M_variablesNames[0] = "v";
    M_variablesNames[1] = "w";
    M_variablesNames[2] = "s";
}

void
FentonKarma4v::setup(GetPot& data, std::string sect)
{
	std::string section = sect + "/FentonKarma4v";
	super::setup(data,sect);
	// From Fenton: Mechanisms of Spiral break up
    // parameter set 1, 3, 4, 5, 6, 8, 9, 10
	// Parameter set 100 for the Right Atria
	// from: A simple model of the right atrium of the human heart with the sinoatrial and atrioventricular nodes included
	//
	M_pulmonary_vein_model = data(section+"/pulmonary_veins",       false);
	if(M_pulmonary_vein_model)
	{
	    std::cout << "* FentonKarma4v: Using pulmonary veins parameters" << std::endl;
	    M_tw1m = 20.0; // ms
	    M_tw2m = 455.0; // ms
	    M_td = 0.17;  // ms
	    M_tsi = 47.5304; // ms
	    M_to = 45.0;  // ms
	    M_uc = 0.25;
	    M_uo = 0.18;
	    M_um = 1.05;
	    M_ucsi = 0.85;
	    M_uso = 0.6;
	    M_k = 5.0;
	    M_aso = 0.025;// ms
	    M_bso = 0.94;// ms
	    M_cso = 0.07;// ms
	}


	double tvp = data(section+"/tvp", -1.0); // ms
	double tv1m = data(section+"/tv1m", -1.0); // ms
	double tv2m = data(section+"/tv2m", -1.0); // ms
	double twp = data(section+"/twp", -1.0); // ms
	double tw1m = data(section+"/tw1m", -1.0); // ms
	double tw2m = data(section+"/tw2m", -1.0); // ms
	double td = data(section+"/td", -1.0);  // ms
	double tsi = data(section+"/tsi", -1.0); // ms
	double to = data(section+"/to", -1.0);  // ms
	double ta = data(section+"/ta", -1.0);  // ms
	double uc = data(section+"/uc", -1.0);
	double uv = data(section+"/uv", -1.0);
	double uw = data(section+"/uw", -1.0);
	double uo = data(section+"/uo", -1.0);
	double um = data(section+"/um", -1.0);
	double ucsi = data(section+"/ucsi", -1.0);
	double uso = data(section+"/uso", -1.0);
	double rsp = data(section+"/rsp", -1.0);
	double rsm = data(section+"/rsm", -1.0);
	double k = data(section+"/k", -1.0);
	double aso = data(section+"/aso", -1.0);// ms
	double bso = data(section+"/bso", -1.0);// ms
	double cso = data(section+"/cso", -1.0);// ms

    M_tvp  = ( tvp  > 0) ? tvp  : M_tvp;
    M_tv1m = ( tv1m > 0) ? tv1m : M_tv1m;
    M_tv2m = ( tv2m > 0) ? tv2m : M_tv2m;
    M_twp  = ( twp  > 0) ? twp  : M_twp;
    M_tw1m = ( tw1m > 0) ? tw1m : M_tw1m;
    M_tw2m = ( tw2m > 0) ? tw2m : M_tw2m;

    M_td  = ( td  > 0) ? td  : M_td;
    M_tsi = ( tsi > 0) ? tsi : M_tsi;
    M_to  = ( to  > 0) ? to  : M_to;
    M_ta  = ( ta  > 0) ? ta  : M_ta;

    M_uc  = ( uc  > 0) ? uc  : M_uc;
    M_uw  = ( uw  > 0) ? uw  : M_uw;
    M_uo  = ( uo  > 0) ? uo  : M_uo;
    M_um  = ( um  > 0) ? um  : M_um;
    M_ucsi= ( ucsi> 0) ? ucsi: M_ucsi;
    M_uso = ( uso > 0) ? uso : M_uso;

    M_rsp = ( rsp > 0) ? rsp : M_rsp;
    M_rsm = ( rsm > 0) ? rsm : M_rsm;
    M_k   = ( k   > 0) ? k   : M_k;
    M_aso = ( aso > 0) ? aso : M_aso;
    M_bso = ( bso > 0) ? bso : M_bso;
    M_cso = ( cso > 0) ? cso : M_cso;
}


void
FentonKarma4v::updateVariables(std::vector<double>& variables, double appliedCurrent, double dt)
{
    double u = variables[0];
    double v = variables[1];
    double w = variables[2];
    double s = variables[3];

    double Humuc = H(u-M_uc);
    double Humuv = H(u-M_uv);
    double Humuw = H(u-M_uw);

    double tvm = M_tv2m * Humuv + M_tv1m * ( 1 - Humuv);
    double twm = M_tw2m * Humuw + M_tw1m * ( 1 - Humuw);
    double rs = M_rsp * Humuc + M_rsm * ( 1 - Humuc);

    double dv = (1 - Humuc) * ( 1 - v ) / tvm - (Humuc) * v / M_tvp;
    double dw = (1 - Humuc) * ( 1 - w ) / twm - (Humuc) * w / M_twp;
    double ds = rs * ( 0.5 * ( 1 + std::tanh( (u-M_ucsi) * M_k ) ) - s );

    // update v
    variables[1] += dt * dv;
    // update w
    variables[2] += dt * dw;
    // update s
    variables[3] += dt * ds;
}

void
FentonKarma4v::updateVariables(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent, double dt)
{
    double u = v_n[0];
    double v = v_n[1];
    double w = v_n[2];
    double s = v_n[3];

    double Humuc = H(u-M_uc);
    double Humuv = H(u-M_uv);
    double Humuw = H(u-M_uw);
    double Hummu = H(M_um-u);

    double tvm = M_tv2m * Humuv + M_tv1m * ( 1 - Humuv);
    double twm = M_tw2m * Humuw + M_tw1m * ( 1 - Humuw);
    double rs = M_rsp * Humuc + M_rsm * ( 1 - Humuc);

    double dv_n = (1 - Humuc) * ( 1 - v ) / tvm - (Humuc) * v / M_tvp;
    double dw_n = (1 - Humuc) * ( 1 - w ) / twm - (Humuc) * w / M_twp;
    double ds_n = rs * ( 0.5 * ( 1 + std::tanh( (u-M_ucsi) * M_k ) ) - s );

    u = v_np1[0];
    v = v_np1[1];
    w = v_np1[2];
    s = v_np1[3];

    Humuc = H(u-M_uc);
    Humuv = H(u-M_uv);
    Humuw = H(u-M_uw);
    Hummu = H(M_um-u);

    tvm = M_tv2m * Humuv + M_tv1m * ( 1 - Humuv);
    twm = M_tw2m * Humuw + M_tw1m * ( 1 - Humuw);
    rs = M_rsp * Humuc + M_rsm * ( 1 - Humuc);

    double dv_np1 = (1 - Humuc) * ( 1 - v ) / tvm - (Humuc) * v / M_tvp;
    double dw_np1 = (1 - Humuc) * ( 1 - w ) / twm - (Humuc) * w / M_twp;
    double ds_np1 = rs * ( 0.5 * ( 1 + std::tanh( (u-M_ucsi) * M_k ) ) - s );


    // update v
    v_np1[1] = v_n[1] + 0.5 * dt * ( dv_n + dv_np1 );
    // update w
    v_np1[2] = v_n[2] + 0.5 * dt * ( dw_n + dw_np1 );
    // update s
    v_np1[3] = v_n[3] + 0.5 * dt * ( ds_n + ds_np1 );
}


void 
FentonKarma4v::updateVariables(std::vector<double>& variables, std::vector<double>& rhs, double appliedCurrent, double dt, bool overwrite)
{
    double u = variables[0];
    double v = variables[1];
    double w = variables[2];
    double s = variables[3];

    double Humuc = H(u-M_uc);
    double Humuv = H(u-M_uv);
    double Humuw = H(u-M_uw);
    double Hummu = H(M_um-u);

    double tvm = M_tv2m * Humuv + M_tv1m * ( 1 - Humuv);
    double twm = M_tw2m * Humuw + M_tw1m * ( 1 - Humuw);
    double rs = M_rsp * Humuc + M_rsm * ( 1 - Humuc);

    rhs[1] = (1 - Humuc) * ( 1 - v ) / tvm - (Humuc) * v / M_tvp;
    rhs[2] = (1 - Humuc) * ( 1 - w ) / twm - (Humuc) * w / M_twp;
    rhs[3] = rs * ( 0.5 * ( 1 + std::tanh( (u-M_ucsi) * M_k ) ) - s );

    if(overwrite)
    {
        variables[1] += dt * rhs[1];
        variables[2] += dt * rhs[2];
        variables[3] += dt * rhs[3];
    }
}



void
FentonKarma4v::updateVariables(double V, std::vector<double>& variables, double dt)
{
    double u = V;
    double v = variables[0];
    double w = variables[1];
    double s = variables[2];

    double Humuc = H(u-M_uc);
    double Humuv = H(u-M_uv);
    double Humuw = H(u-M_uw);
    double Hummu = H(M_um-u);

    double tvm = M_tv2m * Humuv + M_tv1m * ( 1 - Humuv);
    double twm = M_tw2m * Humuw + M_tw1m * ( 1 - Humuw);
    double rs = M_rsp * Humuc + M_rsm * ( 1 - Humuc);

    double dv = (1 - Humuc) * ( 1 - v ) / tvm - (Humuc) * v / M_tvp;
    double dw = (1 - Humuc) * ( 1 - w ) / twm - (Humuc) * w / M_twp;
    double ds = rs * ( 0.5 * ( 1 + std::tanh( (u-M_ucsi) * M_k ) ) - s );

    // update v
    variables[1] += dt * dv;
    // update w
    variables[2] += dt * dw;
    // update s
    variables[3] += dt * ds;

}

double
FentonKarma4v::evaluateIonicCurrent(std::vector<double>& variables, double appliedCurrent, double dt)
{
    double u = variables[0];
    double v = variables[1];
    double w = variables[2];
    double s = variables[3];

    double Humuc = H(u-M_uc);
    double Humuso = H(u-M_uso);

    double Ifi = - v * Humuc * ( u - M_uc) * (M_um - u ) / M_td;

    double tso = M_to;
    if(M_pulmonary_vein_model)
    {
        tso = 2.2 * M_to;
        if( u >= 0.35 ) tso = 1.3 * M_to;
        else if( u >= 0.28) tso = (5.25 - 12.8751 * u) * M_to;
    }
    double Iso = (u - M_uo) * (1 - Humuso) / tso
               + Humuso * M_ta
               + 0.5 * (M_aso - M_ta) * ( 1 + std::tanh( (u - M_bso) / M_cso ) );
    double Isi = - w * s / M_tsi;
    //std::cout << (u - M_uo) << ", " << tso << ", " <<  (1 - Humuso) << std::endl;
    return  Ifi + Iso + Isi;
}


double
FentonKarma4v::evaluateIonicCurrent(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent, double dt)
{
    double u = v_n[0];
    double v = v_n[1];
    double w = v_n[2];
    double s = v_n[3];

    double Humuc = H(u-M_uc);
    double Humuso = H(u-M_uso);

    double Ifi = - v * Humuc * ( u - M_uc) * (M_um - u ) / M_td;

    double tso = M_to;
    if(M_pulmonary_vein_model)
    {
        tso = 2.2 * M_to;
        if( u >= 0.35 ) tso = 1.3 * M_to;
        else if( u >= 0.28) tso = (5.25 - 12.8751 * u) * M_to;
    }
    double Iso = (u - M_uo) * (1 - Humuso) / tso
               + Humuso * M_ta
               + 0.5 * (M_aso - M_ta) * ( 1 + std::tanh( (u - M_bso) / M_cso ) );
    double Isi = - w * s / M_tsi;

    double f_n =   Ifi + Iso + Isi;


    u = v_np1[0];
    v = v_np1[1];
    w = v_np1[2];
    s = v_np1[3];

    Humuc = H(u-M_uc);
    Humuso = H(u-M_uso);

    Ifi = - v * Humuc * ( u - M_uc) * (M_um - u ) / M_td;

    tso = M_to;
    if(M_pulmonary_vein_model)
    {
        tso = 2.2 * M_to;
        if( u >= 0.35 ) tso = 1.3 * M_to;
        else if( u >= 0.28) tso = (5.25 - 12.8751 * u) * M_to;
    }
    Iso = (u - M_uo) * (1 - Humuso) / tso
               + Humuso * M_ta
               + 0.5 * (M_aso - M_ta) * ( 1 + std::tanh( (u - M_bso) / M_cso ) );
    Isi = - w * s / M_tsi;

    double f_np1 =   Ifi + Iso + Isi;

    return 0.5 * (f_n+f_np1);
}

double
FentonKarma4v::evaluateIonicCurrent(double V, std::vector<double>& variables, double appliedCurrent, double dt)
{
    double u = V;
    double v = variables[0];
    double w = variables[1];
    double s = variables[2];

    double Humuc = H(u-M_uc);
    double Humuso = H(u-M_uso);

    double Ifi = - v * Humuc * ( u - M_uc) * (M_um - u ) / M_td;

    double tso = M_to;
    if(M_pulmonary_vein_model)
    {
        tso = 2.2 * M_to;
        if( u >= 0.35 ) tso = 1.3 * M_to;
        else if( u >= 0.28) tso = (5.25 - 12.8751 * u) * M_to;
    }
    double Iso = (u - M_uo) * (1 - Humuso) / tso
               + Humuso * M_ta
               + 0.5 * (M_aso - M_ta) * ( 1 + std::tanh( (u - M_bso) / M_cso ) );
    double Isi = - w * s / M_tsi;

    return  Ifi + Iso + Isi;
}

double
FentonKarma4v::evaluateIonicCurrentTimeDerivative( std::vector<double>& variables,
                                     std::vector<double>& rhs,
                                     double dt,
                                     double h )
{
   return 0;
}

void
FentonKarma4v::initialize(std::vector<double>& variables)
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
FentonKarma4v::initializeSaveData(std::ostream& output)
{
    // time -  0
    output << "time V v w s\n";
}


double
FentonKarma4v::evaluateSAC(double v , double I4f)
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
