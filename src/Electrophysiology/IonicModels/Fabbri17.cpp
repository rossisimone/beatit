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
 * \file Fabbri17.cpp
 *
 * \class Fabbri17
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
 * Created on: Sep 3, 2016
 *
 */

#include "Electrophysiology/IonicModels/Fabbri17.hpp"
#include <cmath>

namespace BeatIt
{

IonicModel* createFabbri17()
{
	return new Fabbri17;
}

Fabbri17::Fabbri17() :
		super(33, 0, "Fabbri17", CellType::MCell)
{
    C = 57e-5; /* cell capacitance (pF) */
    L_cell = 67.0; /* Length of the cell (um) */
    R_cell = 3.9; /* Radius of the cell (um) */

    V_cell = pi * R_cell * R_cell * L_cell; /*cell volume */
    V_sub = 2 * pi * L_sub * L_cell * ( R_cell - 0.5 * L_sub ); /* submembrane space volume */
    V_i = V_ipart * V_cell - V_sub; /*myoplasmic  volume */
    V_jsr = V_jsrpart * V_cell; /*JSR  volume */
    V_nsr = V_nsrpart * V_cell; /*NSR  volume */

	// The name refer to the comments below the lines
	// Nai = 5
    // 1
	M_variablesNames[0] = "Nai"; //
	//y = 0.009508
    // 2
	M_variablesNames[1] = "y"; //
	//m = 0.447724
    // 3
	M_variablesNames[2] = "m";
	// h = 0.003058
    // 4
	M_variablesNames[3] = "h";
	//dL = 0.001921
    // 5
	M_variablesNames[4] = "dL";
	// fL = 0.846702
    // 6
	M_variablesNames[5] = "fL";
	//fCa = 0.844449
    // 7
	M_variablesNames[6] = "fCa";
	// dT = 0.268909
    // 8
	M_variablesNames[7] = "dT";
	//fT = /0.020484
    // 9
	M_variablesNames[8] = "fT";
	//R = 0.9308
    // 10
	M_variablesNames[9] = "R";
	//O = 6.181512e-9
	// 11
	M_variablesNames[10] = "O";
	//I = 4.595622e-10
	// 12
	M_variablesNames[11] = "I";
	//RI = 0.069199
	// 13
	M_variablesNames[12] = "RI";

	//fTC = 0.017929
	// 14
	M_variablesNames[13] = "fTC";
	//fTMC = 0.259947
	// 15
	M_variablesNames[14] = "fTMC";
	//fTMM=0.653777
	// 16
	M_variablesNames[15] = "fTMM";
	//fCMi = 0.217311
	// 17
	M_variablesNames[16] = "fCMi";
	//fCMs = 0.158521
	// 18
	M_variablesNames[17] = "fCMs";
	//fCQ = //0.138975
	// 19
	M_variablesNames[18] = "fCQ";

    // Ca dynamics
	//Cai = 9.15641e-6
    // 20
	M_variablesNames[19] = "Cai";
	//Ca_sub = 6.226104e-5
    // 21
	M_variablesNames[20] = "Ca_sub";
	//Ca_nsr = 0.435148
    // 22
	M_variablesNames[21] = "Ca_nsr";
	//Ca_jsr = 0.409551
    // 23
	M_variablesNames[22] = "Ca_jsr";


    // Ikur
    //r_Kur = 0.011845
    // 24
	M_variablesNames[23] = "r_Kur";
	//s_Kur = 0.845304
    // 25
	M_variablesNames[24] = "s_Kur";


    // i_to
    //q = 0.430836
    // 26
	M_variablesNames[25] = "q";
	// r = 0.014523
    // 27
    M_variablesNames[26] = "r";


    // I_Kr
    // paS = 0.283185
    // 28
    M_variablesNames[27] = "paS";
    // paF = 0.011068
    // 29
    M_variablesNames[28] = "paF";
    // piy = 0.709051
    // 30
    M_variablesNames[29] = "piy";

    // I_Kr
    // paS = 0.1162
    // 31
    M_variablesNames[30] = "n";
    // paF = 0.00277
    // 32
    M_variablesNames[31] = "a";
}

void Fabbri17::initializeSaveData(std::ostream& output)
{
	// time -  0
	output << "time v ";
	for (auto && var_name : M_variablesNames)
	{
		output << var_name << " ";
	}
	output << "\n";
}

void Fabbri17::initialize(std::vector<double>& variables)
{
	// v = -47.787168;       /* Initial Voltage (mv) */
	variables[0] = -47.787168; //
	// Nai = 5(mM)
	variables[1] = 5.0; //
	//y = 0.009508
	variables[2] = 0.009508; //
	//m = 0.447724
	variables[3] = 0.447724;
	// h = 0.003058
	variables[4] = 0.003058;
	//dL=0.001921
	variables[5] = 0.001921;
	// fL = 0.978;
	variables[6] = 0.846702;
	//fCa = 0.000137;
	variables[7] = 0.844449;

    // CaT
	// dT = 0.268909
	variables[8] = 0.268909;
	//fT = 0.020484
	variables[9] = 0.020484;

    //Ca_SR_release
	//R = 0.9308
	variables[10] = 0.9308;
	//O = 6.181512e-9
	variables[11] = 6.181512e-9;
	//I = 4.595622e-10
	variables[12] = 4.595622e-10;
	//RI = 0.069199
	variables[13] = 0.069199;

    // Ca buffering
    // 14
	//fTC = 0.017929
	variables[14] = 0.017929;
	//fTMC = 0.259947
	variables[15] = 0.259947;
	//fTMM=0.653777
	variables[16] = 0.653777;
	//fCMi = 0.217311
	variables[17] = 0.217311;
	//fCMs = 0.158521
	variables[18] = 0.158521;
	//fCQ = 0.138975
	// 19
	variables[19] = 0.138975;


    // Ca dynamics
    // 20
	//Cai = 9.15641e-6
	variables[20] = 9.15641e-6;
	//Ca_sub = 6.226104e-5
	variables[21] = 6.226104e-5;
	//Ca_nsr = 0.435148
	variables[22] = 0.435148;
	//Ca_jsr = 0.409551
	variables[23] = 0.409551;


    // Ikur
	//r_Kur =0.011845
	variables[24] = 0.011845;
	//s_Kur = 0.845304
	variables[25] = 0.845304;


    // i_to
	//q = 0.430836
	variables[26] = 0.430836;
    //r = 0.014523
    variables[27] = 0.014523;

    // I_Kr
    // paS = 0.283185
    variables[28] = 0.283185;
    // paF = 0.011068
    variables[29] = 0.011068;
    // piy = 0.709051
    variables[30] = 0.709051;

    // I_KS
    // n = 0.1162
    variables[31] = 0.1162;
    // I_KACh a gate
    // a = 0.00277
    variables[32] = 0.00277;
}


//! Evaluate total ionic current for the computation of the potential
/*!
 *  \param [in] variables Vector containing the local value of all variables
 *  \param [in] appliedCurrent value of the applied current
 *  \param [in] dt        Timestep
 */
double Fabbri17::evaluateIonicCurrent(std::vector<double>& variables,
		double appliedCurrent, double dt)
{
    // Itot
    i_tot = i_f+i_Kr+i_Ks+i_to+i_NaK+i_NaCa+i_Na+i_CaL+i_CaT+i_KACh+i_Kur;
    //i_tot = i_f+i_Kr+i_Ks+i_to+i_NaK+i_NaCa+i_Na+i_CaL+i_CaT+i_KACh+i_Kur;
//    std::cout << ", i_f: " << i_f;
//    std::cout << ", i_Kr: " << i_Kr;
//    std::cout << ", i_Ks: " << i_Ks;
//    std::cout << ", i_to: " << i_to;
//    std::cout << ", i_NaK: " << i_NaK;
//    std::cout << ", i_NaCa: " << i_NaCa;
//    std::cout << ", i_Na: " << i_Na;
//    std::cout << ", i_CaL: " << i_CaL;
//    std::cout << ", i_CaT: " << i_CaT;
//    std::cout << ", i_KACh: " << i_KACh;
//    std::cout << ", i_Kur: " << i_Kur;
//    std::cout << ", i_tot: " << i_tot << std::endl;
    double dV = i_tot / C;
    //std::cout << "dV: " << dV << ", itot: " << i_tot << ", C: " << C << std::endl;
    return dV;
}

void Fabbri17::updateVariables(std::vector<double>& variables,
		double appliedCurrent, double dt)
{
	this->dt = dt;
	//istim = -appliedCurrent;
	V = variables[0];
	v = variables[0];
	/*  Ion Concentrations */
	Nai = variables[1];
	y = variables[2];

    //Ina
	m = variables[3];
	h = variables[4];
	dL = variables[5];
	fL = variables[6];
	fCa = variables[7];

    // CaT
	dT = variables[8];
	fT = variables[9];

    //Ca_SR_release
	R = variables[10];
	O = variables[11];
	I = variables[12];
	RI = variables[13];

    // Ca buffering
	fTC = variables[14];
	fTMC = variables[15];
	fTMM = variables[16];
	fCMi = variables[17];
	fCMs = variables[18];
	fCQ = variables[19];

    // Ca dynamics
	Cai = variables[20];
	Ca_sub = variables[21];
	Ca_nsr = variables[22];
	Ca_jsr = variables[23];

    // Ikur
	r_Kur = variables[24];
	s_Kur = variables[25];

    // i_to
	q = variables[26];
	r = variables[27];

    // I_Kr
	paS = variables[28];
	paF = variables[29];
	piy = variables[30];

    // I_KS
    n = variables[31];

    // I_KACh a gate
    a = variables[32];

    // UPDATE
    update();

	/*  Ion Concentrations */
	variables[1] = Nai;
	variables[2] = y;
	variables[3] = m;
	variables[4] = h;
	variables[5] = dL;
	variables[6] = fL;
	variables[7] = fCa;

	variables[8] = dT;
	variables[9] = fT;

	variables[10] = R;
	variables[11] = O;
	variables[12] = I;
	variables[13] = RI;

	variables[14] = fTC;
	variables[15] = fTMC;
	variables[16] = fTMM;
	variables[17] = fCMi;
	variables[18] = fCMs;
	variables[19] = fCQ;

	variables[20] = Cai;
	variables[21] = Ca_sub;
	variables[22] = Ca_nsr;
	variables[23] = Ca_jsr;

	variables[24] = r_Kur;
	variables[25] = s_Kur;

	variables[26] = q;
    variables[27] = r;
    variables[28] = paS;
    variables[29] = paF;
    variables[30] = piy;
    variables[31] = n;
    variables[32] = a;

}




/********************************************************/
void
Fabbri17::update()
{
    //def comp Rate_modulation_experiments as
    double  ACh=0.0;//: millimolar {init: 0, pub: out};
    double  Iso_1_uM=0.0;//: dimensionless {init: 0, pub: out};
    // Reversal potentials
     E_Na = RTONF * std::log(Nao/Nai); //reversal potential for Na+
     E_mh = RTONF * std::log( ( Nao+0.12*Ko) / (Nai+0.12*Ki) ); //reversal potential for fast Na+ channel
     E_K  = RTONF * std::log(Ko/Ki); //reversal potential for K+
     E_Ks = RTONF * std::log( ( Ko+0.12*Nao) / (Ki+0.12*Nai) );//reversal potential for slow rectifier K+ channel

    // in paper
    //double E_Ca = RTONF * std::log(Cao/Ca_sub);//reversal potential for Ca2+
    // in code
    double E_Ca = 0.5 * RTONF * std::log(Cao/Ca_sub);//reversal potential for Ca2+

    // Nai_concentration
    //double dNai = (1-Nai_clamp)*-1.0*(i_Na+i_fNa+i_siNa+3.0*i_NaK+3*i_NaCa)/((V_i+V_sub)*F);
    double dNai = -(i_Na+i_fNa+i_siNa+3.0*i_NaK+3*i_NaCa)/((V_i+V_sub)*F);
    Nai += dt * dNai;


    // Funny Current

    // y gating
    y_shift= 0;

    double ACh_shift = 0.0;
    if (ACh > 0 ) ACh_shift = -1-9.898*std::pow(ACh, 0.618)/(std::pow(ACh, 0.618)+0.00122423);

    double Iso_shift = 0.0;
    if (Iso_1_uM > 0 ) Iso_shift = 7.5; //mV

    double tau_y = 1/(0.36*(V+148.8-ACh_shift-Iso_shift)/(exp(0.066*(V+148.8-ACh_shift-Iso_shift))-1)+0.1*(V+87.3-ACh_shift-Iso_shift)/(1-exp(-0.2*(V+87.3-ACh_shift-Iso_shift))))-0.054;

    if ( V < -(80-ACh_shift-Iso_shift-y_shift))
    {
        y_infinity = 0.01329+0.99921/(1+exp((V+97.134-ACh_shift-Iso_shift-y_shift)/8.1752));

    }
    else
    {
        y_infinity = 0.0002501*std::exp(-(V-ACh_shift-Iso_shift-y_shift)/12.861);
    }

    double dy = (y_infinity-y)/tau_y;
    y = y_infinity - (y_infinity - y) * exp(-dt / tau_y);

    // funny current
    double G_f = g_f/(Ko/(Ko+Km_f));
    double G_f_K = G_f/(alpha+1);
    double G_f_Na = alpha*G_f_K;
    double g_f_Na = G_f_Na*Ko/(Ko+Km_f);
    double g_f_K = G_f_K*Ko/(Ko+Km_f);
    double i_fNa = y*g_f_Na*(V-E_Na)*(1-blockade);
    double i_fK = y*g_f_K*(V-E_K)*(1-blockade);
    i_f = i_fNa+i_fK;

    // i_NaK
    double Iso_increase = 1;
    if(Iso_1_uM > 0) Iso_increase = 1.2;

    i_NaK = Iso_increase*i_NaK_max / (1+pow(Km_Kp/Ko, 1.2) )
          / (1+pow(Km_Nap/Nai, 1.3) )
          / (1+exp(-(V-E_Na+110)/20 ) );

    // i_NaCa
    di = 1+Ca_sub/Kci*(1+exp(-Qci*V/RTONF)+Nai/Kcni)+Nai/K1ni*(1+Nai/K2ni*(1+Nai/K3ni));
    doo = 1.0+Cao/Kco*(1+exp(Qco*V/RTONF))+Nao/K1no*(1+Nao/K2no*(1+Nao/K3no));
    double k43 = Nai/(K3ni+Nai);
    double k12 = Ca_sub/Kci*exp(-Qci*V/RTONF)/di;
    double k14 = Nai/K1ni*Nai/K2ni*(1+Nai/K3ni)*exp(Qn*V/(2*RTONF))/di;
    double k41 = exp(-Qn*V/(2*RTONF));
    double k34 = Nao/(K3no+Nao);
    double k21 = Cao/Kco*exp(Qco*V/RTONF)/doo;
    double k23 = Nao/K1no*Nao/K2no*(1+Nao/K3no)*exp(-Qn*V/(2*RTONF))/doo;
    double k32 = exp(Qn*V/(2*RTONF));

    double x1 = k41*k34*(k23+k21)+k21*k32*(k43+k41);
    double x2 = k32*k43*(k14+k12)+k41*k12*(k34+k32);
    double x3 = k14*k43*(k23+k21)+k12*k23*(k43+k41);
    double x4 = k23*k34*(k14+k12)+k14*k21*(k34+k32);

    blockade_NaCa = 0.0;

    i_NaCa = (1-blockade_NaCa)*K_NaCa*(x2*k21-x1*k12)/(x1+x2+x3+x4);

    // i_Na
    i_Na = g_Na*pow(m, 3)*h*(V-E_mh);
    double i_Na_L = g_Na_L*pow(m, 3)*(V-E_mh);
    i_Na = i_Na+i_Na_L;

    // I_Na m gate
    double m_infinity = 1.0/(1+exp(-(V+42.0504)/8.3106));
    double E0_m = V+41;
    double alpha_m = 2000;
    if( abs(E0_m) >= delta_m )
    {
        alpha_m = 200*E0_m/(1-exp(-0.1*E0_m));
    }
    double beta_m = 8000.0*exp(-0.056*(V+66));
    double tau_m = 1.0/(alpha_m+beta_m);
    m = m_infinity - (m_infinity - m) * exp(-dt / tau_m);

    // i_Na_h_gate
    double h_infinity = 1.0/(1+std::exp((V+69.804)/4.4565));
    double alpha_h = 20*exp(-0.125*(V+75));
    double beta_h = 2000.0/(320*exp(-0.1*(V+75))+1.0);
    double tau_h = 1.0/(alpha_h+beta_h);
    h = h_infinity - (h_infinity - h) * exp(-dt / tau_h);

    // i_CaL
    Iso_increase = 1;
    if(Iso_1_uM > 0) Iso_increase = 1.23;

//    std::cout << V << ", " << P_CaL << ", " << RTONF << ", " << Ca_sub << ", " << Cao << ", " << dL << ", " << fL << ", " << fCa << std::endl;
    double i_siCa = 2*P_CaL*(V-0)/(RTONF*(1-exp(-1*(V-0)*2/RTONF)))*(Ca_sub-Cao*exp(-2*(V-0)/RTONF))*dL*fL*fCa;
    double i_siK = 0.000365*P_CaL*(V-0)/(RTONF*(1-exp(-1*(V-0)/RTONF)))*(Ki-Ko*exp(-1*(V-0)/RTONF))*dL*fL*fCa;
    double i_siNa = 0.0000185*P_CaL*(V-0)/(RTONF*(1-exp(-1*(V-0)/RTONF)))*(Nai-Nao*exp(-1*(V-0)/RTONF))*dL*fL*fCa;
    double ACh_block = 0.31*ACh/(ACh+0.00009);
//    std::cout << i_siCa << ", " << i_siK << ", " << i_siNa << ", " << ACh_block << ", " << Iso_increase << std::endl;
    i_CaL = (i_siCa+i_siK+i_siNa)*(1-ACh_block)*1.0*Iso_increase;

    // i_CaL_dL_gate
      double Iso_shift_dL = 0;
      if(Iso_1_uM > 0) Iso_shift_dL = -8;

      double Iso_slope_dL = 0;
      if(Iso_1_uM > 0) Iso_slope_dL = -27;


      double adVm = V;
      if(V == -41.8) V = -41.80001;
      if(V == -6.8) V = -6.80001;

      double bdVm = V;
      if(V == -1.8) V = -1.80001;

      double dL_infinity = 1.0/(1.0+exp(-(V-V_dL-Iso_shift_dL)/(k_dL*(1+Iso_slope_dL/100))));
      double alpha_dL = -0.02839*(adVm+41.8)/(exp(-(adVm+41.8)/2.5)-1)-0.0849*(adVm+6.8)/(exp(-(adVm+6.8)/4.8)-1);
      double beta_dL = 0.01143*(bdVm+1.8)/(exp((bdVm+1.8)/2.5)-1.0);

      double tau_dL = 0.001/(alpha_dL+beta_dL);

      dL = dL_infinity - (dL_infinity - dL) * exp(-dt / tau_dL);

      // i_CaL_fL_gate
      double shift_fL = 0.0;
      double k_fL = 0.0;
      double fL_infinity = 1.0/(1.0+exp((V+37.4+shift_fL)/(5.3+k_fL)));
      double tau_fL = 0.001*(44.3+230.0*exp( -((V+36)/10.0)*((V+36)/10.0) ) );
 //     std::cout << fL_infinity << ", " << tau_fL << std::endl;
      fL = fL_infinity - (fL_infinity - fL) * exp(-dt / tau_fL);

      // i_CaL_fCa_gate
      double fCa_infinity = Km_fCa/(Km_fCa+Ca_sub);
      double tau_fCa = 0.001*fCa_infinity/alpha_fCa;
      fCa = fCa_infinity - (fCa_infinity - fCa) * exp(-dt / tau_fCa);

      // i_CaT
      i_CaT = 2.0*P_CaT*V/(RTONF*(1.0-exp(-1.0*V*2.0/RTONF)))*(Ca_sub-Cao*exp(-2.0*V/RTONF))*dT*fT;

      // i_CaT_dT_gate
      double dT_infinity = 1.0/(1+exp(-(V+38.3)/5.5));
      double tau_dT = 0.001/(1.068*exp((V+38.3)/30)+1.068*exp(-(V+38.3)/30));
      dT = dT_infinity - (dT_infinity - dT) * exp(-dt / tau_dT);

      // i_CaT_fT_gate
      double offset_fT = 0.0;
      double fT_infinity = 1.0/(1+exp((V+58.7)/3.8));
      double tau_fT = 1.0/(16.67*exp(-(V+75)/83.3)+16.67*exp((V+75)/15.38))+offset_fT;
      fT = fT_infinity - (fT_infinity - fT) * exp(-dt / tau_fT);

      // Ca_SR_release
      j_SRCarel = ks*O*(Ca_jsr-Ca_sub);
      double diff = Ca_jsr-Ca_sub;
      double kCaSR = MaxSR-(MaxSR-MinSR)/(1+pow(EC50_SR/Ca_jsr, HSR));
      double koSRCa = koCa/kCaSR;
      double kiSRCa = kiCa*kCaSR;
      double dR = kim*RI-kiSRCa*Ca_sub*R-(koSRCa*Ca_sub*Ca_sub*R-kom*O);
      double dO = koSRCa*Ca_sub*Ca_sub*R-kom*O-(kiSRCa*Ca_sub*O-kim*I);
      double dI = kiSRCa*Ca_sub*O-kim*I-(kom*I-koSRCa*Ca_sub*Ca_sub*RI);
      double dRI = kom*I-koSRCa*Ca_sub*Ca_sub*RI-(kim*RI-kiSRCa*Ca_sub*R);
      R += dt * dR;
      O += dt * dO;
      I += dt * dI;
      RI+= dt * dRI;
      P_tot = R+O+I+RI;

      // Ca_intracellular_fluxes
      double b_up = 0;
      if(Iso_1_uM > 0 ) b_up = -0.25;
      if(ACh > 0 ) b_up =  0.7*ACh/(0.00009+ACh);

    double P_up = P_up_basal*(1-b_up);
    double j_Ca_dif = (Ca_sub-Cai)/tau_dif_Ca;
    double j_up = P_up/(1+exp((-Cai+K_up)/slope_up));
    double j_tr = (Ca_nsr-Ca_jsr)/tau_tr;

    // Ca_buffering
    double delta_fTC = kf_TC*Cai*(1-fTC)-kb_TC*fTC;
    fTC += dt * delta_fTC;

    double delta_fTMC = kf_TMC*Cai*(1-(fTMC+fTMM))-kb_TMC*fTMC;
    fTMC += dt * delta_fTMC;

    double delta_fTMM = kf_TMM*Mgi*(1-(fTMC+fTMM))-kb_TMM*fTMM;
    fTMM += dt * delta_fTMM;

    double delta_fCMi = kf_CM*Cai*(1-fCMi)-kb_CM*fCMi;
    fCMi += dt * delta_fCMi;

    double delta_fCMs = kf_CM*Ca_sub*(1-fCMs)-kb_CM*fCMs;
    fCMs += dt * delta_fCMs;

    double delta_fCQ = kf_CQ*Ca_jsr*(1-fCQ)-kb_CQ*fCQ;
    fCQ += dt * delta_fCQ;

    // Ca_dynamics
    double dCai = 1*(j_Ca_dif*V_sub-j_up*V_nsr)/V_i-(CM_tot*delta_fCMi+TC_tot*delta_fTC+TMC_tot*delta_fTMC);
    double dCa_sub = j_SRCarel*V_jsr/V_sub-((i_siCa+i_CaT-2*i_NaCa)/(2*F*V_sub)+j_Ca_dif+CM_tot*delta_fCMs);
    double dCa_nsr = j_up-j_tr*V_jsr/V_nsr;
    double dCa_jsr = j_tr-(j_SRCarel+CQ_tot*delta_fCQ);

    // i_Kur
    i_Kur = g_Kur*r_Kur*s_Kur*(V-E_K);

    // i_Kur_rKur_gate
    double r_Kur_infinity = 1.0/(1+exp((V+6)/-8.6));
    double tau_r_Kur = 0.009/(1+exp((V+5)/12))+0.0005;
    r_Kur = r_Kur_infinity - (r_Kur_infinity - r_Kur) * exp(-dt / tau_r_Kur);

    // i_Kur_sKur_gate
    double s_Kur_infinity = 1.0/(1.0+exp((V+7.5)/10.0));
    double tau_s_Kur = 0.59/(1+exp((V+60)/10.0))+3.05;
    s_Kur = s_Kur_infinity - (s_Kur_infinity - s_Kur) * exp(-dt / tau_s_Kur);

    // i_to
    i_to = g_to*(V-E_K)*q*r;

    //i_to_q_gate
    double q_infinity = 1.0/(1.0+exp((V+49.0)/13));
    double tau_q = 0.001*0.6*(65.17/(0.57*exp(-0.08*(V+44))+0.065*exp(0.1*(V+45.93)))+10.1);
    q = q_infinity - (q_infinity - q) * exp(-dt / tau_q);

    // i_to_r_gate
    double tau_r = 0.001*0.66*1.4*(15.59/(1.037*exp(0.09*(V+30.61))+0.369*exp(-0.12*(V+23.84)))+2.98);
    double r_infinity = 1.0/(1.0+exp(-(V-19.3)/15));
    r = r_infinity - (r_infinity - r) * exp(-dt / tau_r);

    // i_Kr
    i_Kr = g_Kr*(V-E_K)*(0.9*paF+0.1*paS)*piy;

    // i_Kr_pa_gate
    double alfapaF = 1./(1.+exp(-(V+23.2)/6.6))/(0.84655354/(37.2*exp(V/11.9)+0.96*exp(-V/18.5)));
    double betapaF = 4.*((37.2*exp(V/15.9)+0.96*exp(-V/22.5))/0.84655354-1/(1+exp(-(V+23.2)/10.6))/(0.84655354/(37.2*exp(V/15.9)+0.96*exp(-V/22.5))));
    double pa_infinity = 1/(1+exp(-(V+10.0144)/7.6607));
    double tau_paS = 0.84655354/(4.2*exp(V/17)+0.15*exp(-V/21.6));
    double tau_paF = 1/(30*exp(V/10)+exp(-V/12));
    paS = pa_infinity - (pa_infinity - paS) * exp(-dt / tau_paS);
    paF = pa_infinity - (pa_infinity - paF) * exp(-dt / tau_paF);

    //i_Kr_pi_gate
    double tau_pi = 1/(100*exp(-V/54.645)+656*exp(V/106.157));
    double pi_infinity = 1/(1+exp((V+28.6)/17.1));
    piy = pi_infinity - (pi_infinity - piy) * exp(-dt / tau_pi);

    // i_Ks
    double g_Ks = g_Ks_;
    if (Iso_1_uM > 0) g_Ks = 1.2*g_Ks_;
    i_Ks = g_Ks*(V-E_Ks)*n*n;


    // I_Ks n gate
    Iso_shift = 0.0;
    if(Iso_1_uM > 0) Iso_shift = -14.0;

    double alpha_n = 28.0/(1.0+exp(-(V-40-Iso_shift)/3.0));
    double beta_n = exp(-(V-Iso_shift-5)/25.0);
    double tau_n = 1.0/(alpha_n+beta_n);
    double n_infinity = std::sqrt(1.0/(1.0+exp(-(V+0.6383-Iso_shift)/10.7071)));
    n = n_infinity - (n_infinity - n) * exp(-dt / tau_n);

    // i_KACh
    int ACh_on = 1;
    i_KACh = 0.0;
    if( ACh > 0)
        i_KACh= ACh_on*g_KACh*(V-E_K)*(1+exp((V+20)/20.0))*a;

    // I_KACh a gate
    double alpha_a = (3.5988-0.025641)/(1.0+0.0000012155/pow(ACh, 1.6951))+0.025641;
    double beta_a = 10.0*exp(0.0133*(V+40));
    double tau_a = 1.0/(alpha_a+beta_a);
    double a_infinity = alpha_a/(alpha_a+beta_a);
    a = a_infinity - (a_infinity - a) * exp(-dt / tau_a);


    // Currents:
    // i_tot = i_f+i_Kr+i_Ks+i_to+i_NaK+i_NaCa+i_Na+i_CaL+i_CaT+i_KACh+i_Kur;

}


double
Fabbri17::evaluateIonicCurrentTimeDerivative(std::vector<double>& variables,
        std::vector<double>& old_variables, double dt,
        double h )
{
//    //istim = -appliedCurrent;
//    double qn = old_variables[0];
////    v = variables[0];
////    /*  Ion Concentrations */
//    double dnai = ( variables[1] - old_variables[1] ) / dt;
////    ki = variables[2];
//    double dki = ( variables[2] - old_variables[2] ) / dt;
////    cai = variables[3];
//    double dcai = ( variables[3] - old_variables[3] ) / dt;
////
////    /*  Gate Conditions */
////    m = variables[4];
//    double dm = ( variables[4] - old_variables[4] ) / dt;
////    h = variables[5];
//    double dh = ( variables[5] - old_variables[5] ) / dt;
////    j = variables[6];
//    double dj = ( variables[6] - old_variables[6] ) / dt;
////    d = variables[7];
//    double dd = ( variables[7] - old_variables[7] ) / dt;
////    f = variables[8];
//    double df = ( variables[8] - old_variables[8] ) / dt;
////    xs = variables[9];
//    double dxs = ( variables[9] - old_variables[9] ) / dt;
////    xr = variables[10];
//    double dxr = ( variables[10] - old_variables[10] ) / dt;
////    ato = variables[11];
//    double dato = ( variables[11] - old_variables[11] ) / dt;
////    iito = variables[12];
//    double diito = ( variables[12] - old_variables[12] ) / dt;
////    uakur = variables[13];
//    double duakur = ( variables[13] - old_variables[13] ) / dt;
////    uikur = variables[14];
//    double duikur = ( variables[14] - old_variables[14] ) / dt;
////    fca = variables[15];
//    double dfca = ( variables[15] - old_variables[15] ) / dt;
////    ireljsrol = variables[16];
//    double direljsrol = ( variables[16] - old_variables[16] ) / dt;
////
////    jsr = variables[17];
//    double djsr = ( variables[17] - old_variables[17] ) / dt;
////    nsr = variables[18];
//    double dnsr = ( variables[18] - old_variables[18] ) / dt;
////    trpn = variables[19];
//    double dtrpn = ( variables[19] - old_variables[19] ) / dt;
////    cmdn = variables[20];
//    double dcmdn = ( variables[20] - old_variables[20] ) / dt;
////    csqn = variables[21];
//    double dcsqn = ( variables[21] - old_variables[21] ) / dt;
////    urel = variables[22];
//    double durel = ( variables[22] - old_variables[22] ) / dt;
////    vrel = variables[23];
//    double dvrel = ( variables[23] - old_variables[23] ) / dt;
////    wrel = variables[24];
//    double dwrel = ( variables[24] - old_variables[24] ) / dt;
////    yach = variables[25];
//    double dyach = ( variables[25] - old_variables[25] ) / dt;
////    iky = variables[26];
//    double diky = ( variables[26] - old_variables[26] ) / dt;
//
//    double dItot = 0.0;
//    dItot += comp_d_ina(qn, dnai, dm, dh, dj);
//    dItot += comp_d_ical(qn, dd, df, dfca);
//    dItot += comp_d_ikr(qn, dki, dxr);
//    dItot += comp_d_ikr(qn, dki, dxs);
//    dItot += comp_d_iki(qn, dki);
//    dItot += comp_d_ikach(qn, dki, dyach);
//    dItot += comp_d_ikur(qn, dki, duakur, duikur);
//    dItot += comp_d_ito(qn, dki, dato, diito);
//    dItot += comp_d_inaca(qn, dnai, dcai);
//    dItot += comp_d_inak(qn, dnai);
//    dItot += comp_d_ipca(dcai);
//    dItot += comp_d_icab(qn, dcai);
//    dItot +=  comp_d_inab( qn,  dnai);
//
//    return dItot;
    return 0.0;
}


} /* namespace BeatIt */
