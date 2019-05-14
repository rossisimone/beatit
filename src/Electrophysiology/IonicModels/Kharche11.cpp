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
 * \file Kharche11.cpp
 *
 * \class Kharche11
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

#include "Electrophysiology/IonicModels/Kharche11.hpp"
#include "libmesh/getpot.h"
#include <cmath>
//#include "Electrophysiology/IonicModels/Kharche11LUT.hpp"

namespace BeatIt
{



IonicModel* createKharche11()
{
	return new Kharche11;
}

Kharche11::Kharche11()
 : super(57, 0, "Kharche11", CellType::MCell)
 , set_resting_conditions(false)
 , markov_iks(false)
{
	// The name refer to the comments below the lines
	M_variablesNames[0] = "m";
	//mo=1.405627e-3;
	M_variablesNames[1] = "h";
	//ho= 9.867005e-1;
	M_variablesNames[2] = "j";
	//jo=9.915620e-1;
	M_variablesNames[3] = "d";
	//do=7.175662e-6;
	M_variablesNames[4] = "f";
	//fo=1.000681;
	M_variablesNames[5] = "fcaBj";
	//fcaBjo=2.421991e-2;
	M_variablesNames[6] = "fcaBsl";
	//fcaBslo=1.452605e-2;
	M_variablesNames[7] = "xtof";
	//xtofo=4.051574e-3;
	M_variablesNames[8] = "ytof";
	//ytofo= 9.945511e-1;
	M_variablesNames[9] = "xkr";
	//xkro=8.641386e-3;
	M_variablesNames[10] = "xks";
	//xkso= 5.412034e-3;
	M_variablesNames[11] = "RyRr";
	//RyRro=8.884332e-1;
	M_variablesNames[12] = "RyRo";
	//RyRoo=8.156628e-7;
	M_variablesNames[13] = "RyRi";
	//RyRio=1.024274e-7;
	M_variablesNames[14] = "NaBj";
	//NaBjo=3.539892;
	M_variablesNames[15] = "NaBsl";
	//NaBslo=7.720854e-1;
	M_variablesNames[16] = "TnCL";
	//TnCLo=8.773191e-3;
	M_variablesNames[17] = "TnCHc";
	//TnCHco=1.078283e-1;
	M_variablesNames[18] = "TnCHm";
	//TnCHmo=1.524002e-2;
	M_variablesNames[19] = "CaM";
	//CaMo=2.911916e-4;
	M_variablesNames[20] = "Myoc";
	//Myoco=1.298754e-3;
	M_variablesNames[21] = "Myom";
	//Myomo=1.381982e-1;
	M_variablesNames[22] = "SRB";
	//SRBo=2.143165e-3;
	M_variablesNames[23] = "SLLj";
	//SLLjo=9.566355e-3;
	M_variablesNames[24] = "SLLsl";
	//SLLslo=1.110363e-1;
	M_variablesNames[25] = "SLHj";
	//SLHjo=7.347888e-3;
	M_variablesNames[26] = "SLHsl";
	//SLHslo=7.297378e-2;
	M_variablesNames[27] = "Csqnb";
	//Csqnbo= 1.242988;
	M_variablesNames[28] = "Ca_sr";
	//Ca_sro=0.1e-1; %5.545201e-1;
	M_variablesNames[29] = "Naj";
	//Najo=9.136;%8.80329;
	M_variablesNames[30] = "Nasl";
	//Naslo=9.136;%8.80733;
	M_variablesNames[31] = "Nai";
	//Naio=9.136;%8.80853;
	M_variablesNames[32] = "Ki";
	//Kio=120;
	M_variablesNames[33] = "Caj";
	//Cajo=1.737475e-4;
	M_variablesNames[34] = "Casl";
	//Caslo= 1.031812e-4;
	M_variablesNames[35] = "Cai";
	//Caio=8.597401e-5;
	M_variablesNames[36] = "C1";
	//C1o=0.0015;       % []
	M_variablesNames[37] = "C2";
	//C2o=0.0244;       % []
	M_variablesNames[38] = "C3";
	//C3o=0.1494;       % []
	M_variablesNames[39] = "C4";
	//C4o=0.4071;       % []
	M_variablesNames[40] = "C5";
	//C5o=0.4161;       % []
	M_variablesNames[41] = "C6";
	//C6o=1-(C1o+C2o+C3o+C4o+C5o+C7o+C8o+C9o+C10o+C11o+C12o+C13o+C14o+C15o+O1o+O2o);
	M_variablesNames[42] = "C7";
	//C7o=0.0001;       % []
	M_variablesNames[43] = "C8";
	//C8o=0.0006;       % []
	M_variablesNames[44] = "C9";
	//C9o=0.0008;       % []
	M_variablesNames[45] = "C10";
	//C10o=0;           % []
	M_variablesNames[46] = "C11";
	//C11o=0;           % []
	M_variablesNames[47] = "C12";
	//C12o=0;           % []
	M_variablesNames[48] = "C13";
	//C13o=0;           % []
	M_variablesNames[49] = "C14";
	//C14o=0;           % []
	M_variablesNames[50] = "C15";
	//C15o=0;           % []
	M_variablesNames[51] = "O1";
	//O1o=0;            % []
	M_variablesNames[52] = "rkur";
	//rkuro = 0;
	M_variablesNames[53] = "skur";
	//skuro = 1.0;
	M_variablesNames[54] = "Inal1";
	// INaL1=1
	M_variablesNames[55] = "Inal2";
	// INaL2=0.0
}

void Kharche11::initializeSaveData(std::ostream& output)
{
	// time -  0
	output << "time v ";
	for (auto && var_name : M_variablesNames)
	{
		output << var_name << " ";
	}
	output << "\n";
}

void Kharche11::initialize(std::vector<double>& variables)
{
	if (set_resting_conditions)
	{
		//Vmo=-8.09763e+1;
		variables[0] = -67.13949999999;
		//mo=1.405627e-3;
		variables[1] = 0.0588639;
		//ho= 9.867005e-1;
		variables[2] = 0.1266000;
		//jo=9.915620e-1;
		variables[3] = 0.1266000;
		//do=7.175662e-6;
		variables[4] = 6.19004999999e-05;
		//fo=1.000681;
		variables[5] = 0.9956310;
		//fcaBjo=2.421991e-2;
		variables[6] = 0.0275350;
		//fcaBslo=1.452605e-2;
		variables[7] = 0.0207641;
		// NOTE: xtos and ytos are not used in the original code
		//       I'm skipping thos
		// xtoso=4.051574e-3;
		// variables[8] = 1.2;
		// ytoso=9.945511e-1;
		// variables[9] = 0;
		//xtofo=4.051574e-3;
		variables[8] = 0.002441530;
		//ytofo= 9.945511e-1;
		variables[9] = 0.9102330;
		//xkro=8.641386e-3;
		variables[10] = 0.000010887300;
		//xkso= 5.412034e-3;
		variables[11] = 0.011602300;
		//RyRro=8.884332e-1;
		variables[12] = 0.8179780;
		//RyRoo=8.156628e-7;
		variables[13] = 0.000000477019;
		//RyRio=1.024274e-7;
		variables[14] = 0.0000001061490;
		//NaBjo=3.539892;
		variables[15] = 2.9499700;
		//NaBslo=7.720854e-1;
		variables[16] = 0.64371500;
		//TnCLo=8.773191e-3;
		variables[17] = 0.012254400;
		//TnCHco=1.078283e-1;
		variables[18] = 0.11649400;
		//TnCHmo=1.524002e-2;
		variables[19] = 0.01114040;
		//CaMo=2.911916e-4;
		variables[20] = 0.0004283250;
		//Myoco=1.298754e-3;
		variables[21] = 0.0019062400;
		//Myomo=1.381982e-1;
		variables[22] = 0.1375940;
		//SRBo=2.143165e-3;
		variables[23] = 0.002991050;
		//SLLjo=9.566355e-3;
		variables[24] = 0.008330610;
		//SLLslo=1.110363e-1;
		variables[25] = 0.01372160;
		//SLHjo=7.347888e-3;
		variables[26] = 0.079161100;
		//SLHslo=7.297378e-2;
		variables[27] = 0.144150;
		//Csqnbo= 1.242988;
		variables[28] = 0.82474700;
		//Ca_sro=0.1e-1; %5.545201e-1;
		variables[29] = 0.30197700;
		//Najo=9.136;%8.80329;
		variables[30] = 6.397640;
		//Naslo=9.136;%8.80733;
		variables[31] = 6.396940;
		//Naio=9.136;%8.80853;
		variables[32] = 6.396940;
		//Kio=120;
		variables[33] = 120;
		//Cajo=1.737475e-4;
		variables[34] = 0.0001982020; // Assume Cai diastolic
		//Caslo= 1.031812e-4;
		variables[35] = 0.000148430;
		//Caio=8.597401e-5;
		variables[36] = 0.00012719800;
		// NOTE: rtos is available in the provided matlab code but it's not used anywhere
		//       For the time being I will comment it
		// rtoso=0.9946;
		// variables[39] = 0;
		//C1o=0.0015;       % []
		double sum = 0.0;
		variables[37] = 0.0015;
		sum += variables[37];
		//C2o=0.0244;       % []
		variables[38] = 0.0244;
		sum += variables[38];
		//C3o=0.1494;       % []
		variables[39] = 0.1494;
		sum += variables[39];
		//C4o=0.4071;       % []
		variables[40] = 0.4071;
		sum += variables[40];
		//C5o=0.4161;       % []
		variables[41] = 0.4161;
		sum += variables[41];
		// NOTE: C6 is put later so we skip one index: index 42
		//C7o=0.0001;       % []
		variables[43] = 0.0001;
		sum += variables[43];
		//C8o=0.0006;       % []
		variables[44] = 0.0006;
		sum += variables[44];
		//C9o=0.0008;       % []
		variables[45] = 0.0008;
		sum += variables[45];
		//C10o=0;           % []
		variables[46] = 0.0;
		//C11o=0;           % []
		variables[47] = 0.0;
		//C12o=0;           % []
		variables[48] = 0.0;
		//C13o=0;           % []
		variables[49] = 0.0;
		//C14o=0;           % []
		variables[50] = 0.0;
		//C15o=0;           % []
		variables[51] = 0.0;
		//O1o=0;            % []
		variables[52] = 0.0;
		//O2o=0;            % []
		// NOTE: C6 is set here we fill here the missing index 42
		double O2o = 0.0;
		//C6o=1-(C1o+C2o+C3o+C4o+C5o+C7o+C8o+C9o+C10o+C11o+C12o+C13o+C14o+C15o+O1o+O2o);
		variables[42] = 1 - sum;
		//initial conditions for IKur
		//rkuro = 0;
		variables[53] = 8.16842e-04;
		//skuro = 1.0;
		variables[54] = 0.9974370;
		// For INaL
		// y70=[1; 0; 0];
		// NOTE: there are no names here ...
		//        I'll call them INal1 and INal2
		//        INal3 is not used so I will not use it
		// INaL1=1
		variables[55] = 0.02720480;
		// INaL1=0.0
		variables[56] = 0.0196170;
	}
	else
	{
		//Vmo=-8.09763e+1;
		variables[0] = -80.9763;
		//mo=1.405627e-3;
		variables[1] = 1.405627e-3;
		//ho= 9.867005e-1;
		variables[2] = 9.867005e-1;
		//jo=9.915620e-1;
		variables[3] = 9.915620e-1;
		//do=7.175662e-6;
		variables[4] = 7.175662e-6;
		//fo=1.000681;
		variables[5] = 1.000681;
		//fcaBjo=2.421991e-2;
		variables[6] = 2.421991e-2;
		//fcaBslo=1.452605e-2;
		variables[7] = 1.452605e-2;
		// NOTE: xtos and ytos are not used in the original code
		//       I'm skipping thos
		// xtoso=4.051574e-3;
		// variables[8] = 1.2;
		// ytoso=9.945511e-1;
		// variables[9] = 0;
		//xtofo=4.051574e-3;
		variables[8] = 4.051574e-3;
		//ytofo= 9.945511e-1;
		variables[9] = 9.945511e-1;
		//xkro=8.641386e-3;
		variables[10] = 8.641386e-3;
		//xkso= 5.412034e-3;
		variables[11] = 5.412034e-3;
		//RyRro=8.884332e-1;
		variables[12] = 8.884332e-1;
		//RyRoo=8.156628e-7;
		variables[13] = 8.156628e-7;
		//RyRio=1.024274e-7;
		variables[14] = 1.024274e-7;
		//NaBjo=3.539892;
		variables[15] = 3.539892;
		//NaBslo=7.720854e-1;
		variables[16] = 7.720854e-1;
		//TnCLo=8.773191e-3;
		variables[17] = 8.773191e-3;
		//TnCHco=1.078283e-1;
		variables[18] = 1.078283e-1;
		//TnCHmo=1.524002e-2;
		variables[19] = 1.524002e-2;
		//CaMo=2.911916e-4;
		variables[20] = 2.911916e-4;
		//Myoco=1.298754e-3;
		variables[21] = 1.298754e-3;
		//Myomo=1.381982e-1;
		variables[22] = 1.381982e-1;
		//SRBo=2.143165e-3;
		variables[23] = 2.143165e-3;
		//SLLjo=9.566355e-3;
		variables[24] = 9.566355e-3;
		//SLLslo=1.110363e-1;
		variables[25] = 1.110363e-1;
		//SLHjo=7.347888e-3;
		variables[26] = 7.347888e-3;
		//SLHslo=7.297378e-2;
		variables[27] = 7.297378e-2;
		//Csqnbo= 1.242988;
		variables[28] = 1.242988;
		//Ca_sro=0.1e-1; %5.545201e-1;
		variables[29] = 0.1e-1;
		//Najo=9.136;%8.80329;
		variables[30] = 9.136;
		//Naslo=9.136;%8.80733;
		variables[31] = 9.136;
		//Naio=9.136;%8.80853;
		variables[32] = 9.136;
		//Kio=120;
		variables[33] = 120;
		//Cajo=1.737475e-4;
		variables[34] = 1.737475e-4; // Assume Cai diastolic
		//Caslo= 1.031812e-4;
		variables[35] = 1.031812e-4;
		//Caio=8.597401e-5;
		variables[36] = 8.597401e-5;
		// NOTE: rtos is available in the provided matlab code but it's not used anywhere
		//       For the time being I will comment it
		// rtoso=0.9946;
		// variables[39] = 0;
		//C1o=0.0015;       % []
		double sum = 0.0;
		variables[37] = 0.0015;
		sum += variables[37];
		//C2o=0.0244;       % []
		variables[38] = 0.0244;
		sum += variables[38];
		//C3o=0.1494;       % []
		variables[39] = 0.1494;
		sum += variables[39];
		//C4o=0.4071;       % []
		variables[40] = 0.4071;
		sum += variables[40];
		//C5o=0.4161;       % []
		variables[41] = 0.4161;
		sum += variables[41];
		// NOTE: C6 is put later so we skip one index: index 42
		//C7o=0.0001;       % []
		variables[43] = 0.0001;
		sum += variables[43];
		//C8o=0.0006;       % []
		variables[44] = 0.0006;
		sum += variables[44];
		//C9o=0.0008;       % []
		variables[45] = 0.0008;
		sum += variables[45];
		//C10o=0;           % []
		variables[46] = 0.0;
		//C11o=0;           % []
		variables[47] = 0.0;
		//C12o=0;           % []
		variables[48] = 0.0;
		//C13o=0;           % []
		variables[49] = 0.0;
		//C14o=0;           % []
		variables[50] = 0.0;
		//C15o=0;           % []
		variables[51] = 0.0;
		//O1o=0;            % []
		variables[52] = 0.0;
		//O2o=0;            % []
		// NOTE: C6 is set here we fill here the missing index 42
		double O2o = 0.0;
		//C6o=1-(C1o+C2o+C3o+C4o+C5o+C7o+C8o+C9o+C10o+C11o+C12o+C13o+C14o+C15o+O1o+O2o);
		variables[42] = 1 - sum;
		//initial conditions for IKur
		//rkuro = 0;
		variables[53] = 0.0;
		//skuro = 1.0;
		variables[54] = 1.0;
		// For INaL
		// y70=[1; 0; 0];
		// NOTE: there are no names here ...
		//        I'll call them INal1 and INal2
		//        INal3 is not used so I will not use it
		// INaL1=1
		variables[55] = 1.0;
		// INaL1=0.0
		variables[56] = 0.0;
	}
}

void Kharche11::setup(GetPot& data, std::string sect)
{
	std::string section = sect + "/Kharche11";
	set_resting_conditions = data(section + "/resting_values", false);
	std::cout << "Setting resting values: "<< set_resting_conditions << ", " << section << std::endl;
}


//! Evaluate total ionic current for the computation of the potential
/*!
 *  \param [in] variables Vector containing the local value of all variables
 *  \param [in] appliedCurrent value of the applied current
 *  \param [in] dt        Timestep
 */
double Kharche11::evaluateIonicCurrent(std::vector<double>& variables,
		double appliedCurrent, double dt)
{
	// For compatibility  with the original code where the applied stimulus in opposite
    // Do not include applied current
	return I_tot;

}

//Nernst Potentials
void Kharche11::nernst_potentials(const std::vector<double>& variables)
{
	//Najo=9.136;%8.80329;
	double Naj = variables[30];
	//Naslo=9.136;%8.80733;
	double Nasl = variables[31];
	//Kio=120;
	double Ki = variables[33];
	//Cajo=1.737475e-4;
	double Caj = variables[34];
	//Caslo= 1.031812e-4;
	double Casl = variables[35];
	ena_junc = (1. / FoRT) * std::log(Nao / Naj);     // [mV]
	ena_sl = (1. / FoRT) * std::log(Nao / Nasl);       // [mV]
	ek = (1. / FoRT) * std::log(Ko / Ki);            // [mV]
	eca_junc = (1. / FoRT / 2) * std::log(Cao / Caj);   // [mV]
	eca_sl = (1. / FoRT / 2) * std::log(Cao / Casl);     // [mV]
	ecl = (1. / FoRT) * std::log(Cli / Clo);            // [mV]
}

double Kharche11::evaluateIonicCurrentTimeDerivative(std::vector<double>& variables,
		std::vector<double>& old_variables, double dt, double dx)
{
	// Membrane Currents
	double v = variables[0];
	double Q = old_variables[0];

	double m = variables[1];
	double h = variables[2];
	double j = variables[3];
	double d = variables[4];
	double f = variables[5];
	double fcaBj = variables[6];
	double fcaBsl = variables[7];
	double xtof = variables[8];
	double ytof = variables[9];
	double xkr = variables[10];
	double xks = variables[11];
	double RyRr = variables[12];
	double RyRo = variables[13];
	double RyRi = variables[14];
	double NaBj = variables[15];
	double NaBsl = variables[16];
	double TnCL = variables[17];
	double TnCHc = variables[18];
	double TnCHm = variables[19];
	double CaM = variables[20];
	double Myoc = variables[21];
	double Myom = variables[22];
	double SRB = variables[23];
	double SLLj = variables[24];
	double SLLsl = variables[25];
	double SLHj = variables[26];
	double SLHsl = variables[27];
	double Csqnb = variables[28];
	double Ca_sr = variables[29];
	double Naj = variables[30];
	double Nasl = variables[31];
	double Nai = variables[32];
	double Ki = variables[33];
	double Caj = variables[34];
	double Casl = variables[35];
	double Cai = variables[36];
	double C1 = variables[37];
	double C2 = variables[38];
	double C3 = variables[39];
	double C4 = variables[40];
	double C5 = variables[41];
	double C6 = variables[42];
	double C7 = variables[43];
	double C8 = variables[44];
	double C9 = variables[45];
	double C10 = variables[46];
	double C11 = variables[47];
	double C12 = variables[48];
	double C13 = variables[49];
	double C14 = variables[50];
	double C15 = variables[51];
	double O1 = variables[52];
	double rkur = variables[53];
	double skur = variables[54];
	double Inal1 = variables[55];
	double Inal2 = variables[56];

	double dm = (variables[1] - old_variables[1]) / dt;
	double dh = (variables[2] - old_variables[2]) / dt;
	double dj = (variables[3] - old_variables[3]) / dt;
	double dd = (variables[4] - old_variables[4]) / dt;
	double df = (variables[5] - old_variables[5]) / dt;
	double dfcaBj = (variables[6] - old_variables[6]) / dt;
	double dfcaBsl = (variables[7] - old_variables[7]) / dt;
	double dxtof = (variables[8] - old_variables[8]) / dt;
	double dytof = (variables[9] - old_variables[9]) / dt;
	double dxkr = (variables[10] - old_variables[10]) / dt;
	double dxks = (variables[11] - old_variables[11]) / dt;
	double dRyRr = (variables[12] - old_variables[12]) / dt;
	double dRyRo = (variables[13] - old_variables[13]) / dt;
	double dRyRi = (variables[14] - old_variables[14]) / dt;
	double dNaBj = (variables[15] - old_variables[15]) / dt;
	double dNaBsl = (variables[16] - old_variables[16]) / dt;
	double dTnCL = (variables[17] - old_variables[17]) / dt;
	double dTnCHc = (variables[18] - old_variables[18]) / dt;
	double dTnCHm = (variables[19] - old_variables[19]) / dt;
	double dCaM = (variables[20] - old_variables[20]) / dt;
	double dMyoc = (variables[21] - old_variables[21]) / dt;
	double dMyom = (variables[22] - old_variables[22]) / dt;
	double dSRB = (variables[23] - old_variables[23]) / dt;
	double dSLLj = (variables[24] - old_variables[24]) / dt;
	double dSLLsl = (variables[25] - old_variables[25]) / dt;
	double dSLHj = (variables[26] - old_variables[26]) / dt;
	double dSLHsl = (variables[27] - old_variables[27]) / dt;
	double dCsqnb = (variables[28] - old_variables[28]) / dt;
	double dCa_sr = (variables[29] - old_variables[29]) / dt;
	double dNaj = (variables[30] - old_variables[30]) / dt;
	double dNasl = (variables[31] - old_variables[31]) / dt;
	double dNai = (variables[32] - old_variables[32]) / dt;
	double dKi = (variables[33] - old_variables[33]) / dt;
	double dCaj = (variables[34] - old_variables[34]) / dt;
	double dCasl = (variables[35] - old_variables[35]) / dt;
	double dCai = (variables[36] - old_variables[36]) / dt;
	double dC1 = (variables[37] - old_variables[37]) / dt;
	double dC2 = (variables[38] - old_variables[38]) / dt;
	double dC3 = (variables[39] - old_variables[39]) / dt;
	double dC4 = (variables[40] - old_variables[40]) / dt;
	double dC5 = (variables[41] - old_variables[41]) / dt;
	double dC6 = (variables[42] - old_variables[42]) / dt;
	double dC7 = (variables[43] - old_variables[43]) / dt;
	double dC8 = (variables[44] - old_variables[44]) / dt;
	double dC9 = (variables[45] - old_variables[45]) / dt;
	double dC10 = (variables[46] - old_variables[46]) / dt;
	double dC11 = (variables[47] - old_variables[47]) / dt;
	double dC12 = (variables[48] - old_variables[48]) / dt;
	double dC13 = (variables[49] - old_variables[49]) / dt;
	double dC14 = (variables[50] - old_variables[50]) / dt;
	double dC15 = (variables[51] - old_variables[51]) / dt;
	double dO1 = (variables[52] - old_variables[52]) / dt;
	double drkur = (variables[53] - old_variables[53]) / dt;
	double dskur = (variables[54] - old_variables[54]) / dt;
	double dInal1 = (variables[55] - old_variables[55]) / dt;
	double dInal2 = (variables[56] - old_variables[56]) / dt;

	double aux = 0.0;
	// INa
	double dI_Na_junc = Fjunc * GNa * m * m * m * h * j * Q
			+ Fjunc * GNa * 3 * m * m * dm * h * j * (v - ena_junc)
			+ Fjunc * GNa * m * m * m * dh * j * (v - ena_junc)
			+ Fjunc * GNa * m * m * m * h * dj * (v - ena_junc);

	double dI_Na_sl = Fsl * GNa * m * m * m * h * j * Q
			+ Fsl * GNa * 3 * m * m * dm * h * j * (v - ena_sl)
			+ Fsl * GNa * m * m * m * dh * j * (v - ena_sl)
			+ Fsl * GNa * m * m * m * h * dj * (v - ena_sl);
	double dI_Na = dI_Na_junc + dI_Na_sl;

	double dI_NaL_junc = Fjunc * GNaL * Inal1 * Inal1 * Inal1 * Inal2 * Q
			+ Fjunc * GNaL * 3 * Inal1 * Inal1 * dInal1 * Inal2 * (v - ena_junc)
			+ Fjunc * GNaL * Inal1 * Inal1 * Inal1 * dInal2 * (v - ena_junc);

	double dI_NaL_sl = Fsl * GNaL * Inal1 * Inal1 * Inal1 * Inal2 * Q
			+ Fsl * GNaL * 3 * Inal1 * Inal1 * dInal1 * Inal2 * (v - ena_sl)
			+ Fsl * GNaL * Inal1 * Inal1 * Inal1 * dInal2 * (v - ena_sl);
	double dI_NaL = dI_NaL_junc + dI_NaL_sl;

	// I_nabk: Na Background Current
	double dI_nabk_junc = Fjunc * GNaB * Q;
	double dI_nabk_sl = Fsl * GNaB * Q;
	double dI_nabk = dI_nabk_junc + dI_nabk_sl;

	// I_nak: Na/K Pump Current
	double sigma = (std::exp(Nao / 67.3) - 1) / 7.0;
	double fnak = 1.0
			/ (1.0 + 0.1245 * std::exp(-0.1 * v * FoRT)
					+ 0.0365 * sigma * std::exp(-v * FoRT));
	double dfnak = -fnak * fnak
			* (-0.1 * FoRT * 0.1245 * std::exp(-0.1 * v * FoRT)
					- FoRT * 0.0365 * sigma * std::exp(-v * FoRT));
	aux = KmNaip / Naj;
	double daux_Naj = -KmNaip / Naj / Naj * dNaj;

	double dI_nak_junc = Fjunc * IbarNaK * dfnak * Ko
			/ (1.0 + aux * aux * aux * aux) / (Ko + KmKo) * Q
			+ Fjunc * IbarNaK * fnak * Ko / (1.0 + aux * aux * aux * aux)
					/ (1.0 + aux * aux * aux * aux) / (Ko + KmKo) * 4 * aux
					* aux * aux * daux_Naj;

	aux = KmNaip / Nasl;
	double daux_Nasl = -KmNaip / Nasl / Nasl * dNasl;
	double dI_nak_sl = Fsl * IbarNaK * dfnak * Ko
			/ (1.0 + aux * aux * aux * aux) / (Ko + KmKo) * Q
			+ Fsl * IbarNaK * fnak * Ko / (1.0 + aux * aux * aux * aux)
					/ (1.0 + aux * aux * aux * aux) / (Ko + KmKo) * 4 * aux
					* aux * aux * daux_Nasl;

	double dI_nak = dI_nak_junc + dI_nak_sl;

	// I_kr: Rapidly Activating K Current
	double gkr = 0.035 * std::sqrt(Ko / 5.4);
	double rkr = 1.0 / (1.0 + std::exp((v + 74.0) / 24.0));
	double drkr = -rkr * rkr * std::exp((v + 74.0) / 24.0) / 24.0;
	double dI_kr = (gkr * xkr * drkr * (v - ek) + gkr * xkr * rkr) * Q;

	// I_ks: Slowly Activating K Current

	double eks = (1 / FoRT) * std::log((Ko + pNaK * Nao) / (Ki + pNaK * Nai));
	double deks = (1 / FoRT) * (-(dKi + pNaK * dNai)) / (Ki + pNaK * Nai);

	double dI_ks = 0.0;
	if (!markov_iks)
	{
		double gks_junc = 1.0 * (1.0 + 1.0 * AF + 2.0 * ISO) * 0.0035 * 1.0;
		double gks_sl = 1.0 * (1.0 + 1.0 * AF + 2.0 * ISO) * 0.0035 * 1.0; // Fra

		double dI_ks_junc = Fjunc * gks_junc * xks * xks * (Q - deks)
				+ Fjunc * gks_junc * 2 * xks * dxks * (v - eks);
		double dI_ks_sl = Fsl * gks_sl * xks * xks * (Q - deks)
				+ Fsl * gks_sl * 2 * xks * dxks * (v - eks);
		dI_ks = dI_ks_junc + dI_ks_sl;
	}
	else
	{
		double gks_junc = 1 * 0.0065;
		double gks_sl = 1 * 0.0065; //FRA
		double alpha = 3.98e-4 * std::exp(3.61e-1 * v * FoRT);
		double beta = 5.74e-5 * std::exp(-9.23e-2 * v * FoRT);
		double gamma = 3.41e-3 * std::exp(8.68e-1 * v * FoRT);
		double delta = 1.2e-3 * std::exp(-3.3e-1 * v * FoRT);
		double teta = 6.47e-3;
		double eta = 1.25e-2 * std::exp(-4.81e-1 * v * FoRT);
		double psi = 6.33e-3 * std::exp(1.27 * v * FoRT);
		double omega = 4.91e-3 * std::exp(-6.79e-1 * v * FoRT);

		double O2 = 1
				- (C1 + C2 + C3 + C4 + C5 + C6 + C8 + C7 + C9 + C10 + C11 + C12
						+ C13 + C14 + C15 + O1);
		double dO2 = -(dC1 + dC2 + dC3 + dC4 + dC5 + dC6 + dC8 + dC7 + dC9
				+ dC10 + dC11 + dC12 + dC13 + dC14 + dC15 + dO1);

		double dI_ks_junc = Fjunc * gks_junc * (O1 + O2) * (Q - deks)
				+ Fjunc * gks_junc * (dO1 + dO2) * (v - eks);
		double dI_ks_sl = Fsl * gks_sl * (O1 + O2) * (Q - deks)
				+ Fsl * gks_sl * (dO1 + dO2) * (v - eks);
		dI_ks = dI_ks_junc + dI_ks_sl;
	}

	// I_kp: Plateau K current
	double kp_kp = 1.0 / (1.0 + std::exp(7.488 - v / 5.98));
	double dkp_kp = -kp_kp * kp_kp * (std::exp(7.488 - v / 5.98) / (-5.98));
	double dI_kp_junc = Fjunc * gkp * dkp_kp * (v - ek) + Fjunc * gkp * kp_kp;
	double dI_kp_sl = Fsl * gkp * dkp_kp * (v - ek) + Fsl * gkp * kp_kp;
	double dI_kp = (dI_kp_junc + dI_kp_sl) * Q;

	// I_to: Transient Outward K Current (slow and fast components)

	//ytof += dt *( ytoss - ytof )  / tauytof;
	double dI_tof = GtoFast * xtof * ytof * Q
			+ GtoFast * dxtof * ytof * (v - ek)
			+ GtoFast * xtof * dytof * (v - ek);
	double dI_to = dI_tof;

	// I_kur: Ultra rapid delayed rectifier Outward K Current
	// Equation for IKur; from Maleckar et al. 2009 - EG
	// atrium
	// equations for activation;
	double RA = 0; // Right Atrium
	double Gkur = 1 * (1.0 - 0.5 * AF) * (1 + 2 * ISO) * 0.045 * (1 + 0.2 * RA); // nS/pF maleckar 0.045

	double dI_kur = Gkur * rkur * skur * Q + Gkur * drkur * skur * (v - ek)
			+ Gkur * rkur * dskur * (v - ek);

	///
	// I_ki: Time-Independent K Current
	double aki = 1.02 / (1 + std::exp(0.2385 * (v - ek - 59.215)));
	double daki = -aki * aki / 1.02
			* (0.2385 * std::exp(0.2385 * (v - ek - 59.215)));
	double bki = (0.49124 * std::exp(0.08032 * (v + 5.476 - ek))
			+ std::exp(0.06175 * (v - ek - 594.31)))
			/ (1.0 + std::exp(-0.5143 * (v - ek + 4.753)));
	double dbki = (0.06175 * std::exp(0.06175 * v - 0.06175 * ek - 36.6986425)
			+ 0.0394563968 * std::exp(0.08032 * v - 0.08032 * ek + 0.43983232))
			/ (std::exp(0.5143 * ek - 0.5143 * v - 2.4444679) + 1.0)
			+ (0.5143 * std::exp(0.5143 * ek - 0.5143 * v - 2.4444679)
					* (std::exp(0.06175 * v - 0.06175 * ek - 36.6986425)
							+ 0.49124
									* std::exp(
											0.08032 * v - 0.08032 * ek
													+ 0.43983232)))
					/ (std::exp(0.5143 * ek - 0.5143 * v - 2.4444679) + 1.0)
					/ (std::exp(0.5143 * ek - 0.5143 * v - 2.4444679) + 1.0);
	double kiss = aki / (aki + bki);
	double dkiss = daki / (aki + bki)
			- aki / (aki + bki) / (aki + bki) * (daki + dbki);

	// I_ki =1* 0.35*sqrt(Ko/5.4)*kiss*(v-ek);
	// SVP 11/11/09
	// multiplieD IK1 by 0.15 to scale it to single cell isolated atrial cell
	// resting potential
	double dI_ki = (1.0 + 1.0 * AF) * 0.0525 * std::sqrt(Ko / 5.4) * dkiss
			* (v - ek) + (1.0 + 1.0 * AF) * 0.0525 * std::sqrt(Ko / 5.4) * kiss;
	dI_ki *= Q;

	// I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
	double dI_ClCa_junc = Fjunc * GClCa / (1 + KdClCa / Caj) * Q
			+ Fjunc * GClCa / (1 + KdClCa / Caj) / (1 + KdClCa / Caj)
					* (-KdClCa / Caj / Caj * dCaj) * (v - ecl);
	double dI_ClCa_sl = Fsl * GClCa / (1 + KdClCa / Casl) * Q
			+ Fsl * GClCa / (1 + KdClCa / Casl) / (1 + KdClCa / Casl)
					* (-KdClCa / Casl / Casl * dCasl) * (v - ecl);
	double dI_ClCa = dI_ClCa_junc + dI_ClCa_sl;
	double dI_Clbk = GClB * Q;

	double GClCFTR = 0; // 4.9e-3*ISO;     // [mS/uF]
	double dI_ClCFTR = GClCFTR;

	// I_Ca: L-type Calcium Current
	double dss = 1.0 / (1.0 + std::exp(-(v + 3 * ISO + 9) / 6)); //in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
	double taud = 1.0 * dss * (1 - std::exp(-(v + 3 * ISO + 9) / 6))
			/ (0.035 * (v + 3 * ISO + 9));
	double fss = 1.0 / (1.0 + std::exp((v + 3 * ISO + 30) / 7))
			+ 0.2 / (1 + std::exp((50 - v - 3 * ISO) / 20)); // in Maleckar v1/2=-27.4 S=7.1 (mV); Courtemanche v1/2=-28 S=6.9 (mV)
	double tauf = 1.0
			/ (0.0197
					* std::exp(
							-(0.0337 * (v + 3 * ISO + 25))
									* (0.0337 * (v + 3 * ISO + 25))) + 0.02);

	// ////////////////////////////////////////////////////////
	// ////////////////////////////////////////////////////////
	double fcaCaMSL = 0.1 / (1.0 + (0.01 / Casl));
	double dfcaCaMSL = -0.1 / (1.0 + (0.01 / Casl)) / (1.0 + (0.01 / Casl))
			* (-0.01 / Casl / Casl * dCasl);
	double fcaCaj = 0.1 / (1.0 + (0.01 / Caj));
	double dfcaCaj = -0.1 / (1.0 + (0.01 / Caj)) / (1.0 + (0.01 / Caj))
			* (-0.01 / Caj / Caj * dCaj);
	fcaCaMSL = 0.0;
	fcaCaj = 0.0;
	double ibarca_j = pCa * 4 * (v * Frdy * FoRT)
			* (0.341 * Caj * std::exp(2 * v * FoRT) - 0.341 * Cao)
			/ (std::exp(2 * v * FoRT) - 1.0);
	double dibarca_j = (2.728 * Caj * FoRT * FoRT * Frdy * pCa * v
			* exp(2.0 * FoRT * v)) / (std::exp(2.0 * FoRT * v) - 1.0)
			- (4.0 * FoRT * Frdy * pCa
					* (0.341 * Cao - 0.341 * Caj * std::exp(2.0 * FoRT * v)))
					/ (std::exp(2.0 * FoRT * v) - 1.0)
			+ (8.0 * FoRT * FoRT * Frdy * pCa * v * std::exp(2.0 * FoRT * v)
					* (0.341 * Cao - 0.341 * Caj * std::exp(2.0 * FoRT * v)))
					/ (std::exp(2.0 * FoRT * v) - 1.0)
					/ (std::exp(2.0 * FoRT * v) - 1.0);
	dibarca_j *= Q;
	dibarca_j += pCa * 4 * (v * Frdy * FoRT)
			* (0.341 * dCaj * std::exp(2 * v * FoRT))
			/ (std::exp(2 * v * FoRT) - 1.0);

	double ibarca_sl = pCa * 4 * (v * Frdy * FoRT)
			* (0.341 * Casl * std::exp(2 * v * FoRT) - 0.341 * Cao)
			/ (std::exp(2 * v * FoRT) - 1.0);
	double dibarca_sl = (2.728 * Casl * FoRT * FoRT * Frdy * pCa * v
			* std::exp(2.0 * FoRT * v)) / (exp(2.0 * FoRT * v) - 1.0)
			- (4.0 * FoRT * Frdy * pCa
					* (0.341 * Cao - 0.341 * Casl * std::exp(2.0 * FoRT * v)))
					/ (std::exp(2.0 * FoRT * v) - 1.0)
			+ (8.0 * FoRT * FoRT * Frdy * pCa * v * std::exp(2.0 * FoRT * v)
					* (0.341 * Cao - 0.341 * Casl * std::exp(2.0 * FoRT * v)))
					/ (std::exp(2.0 * FoRT * v) - 1.0)
					/ (std::exp(2.0 * FoRT * v) - 1.0);
	dibarca_sl *= Q;
	dibarca_sl += pCa * 4 * (v * Frdy * FoRT)
			* (0.341 * dCasl * std::exp(2 * v * FoRT))
			/ (std::exp(2 * v * FoRT) - 1.0);

	double ibark = pK * (v * Frdy * FoRT)
			* (0.75 * Ki * std::exp(v * FoRT) - 0.75 * Ko)
			/ (std::exp(v * FoRT) - 1.0);
	double dibark = (FoRT * FoRT * Frdy * pK * v * std::exp(FoRT * v)
			* (0.75 * Ko - 0.75 * Ki * std::exp(FoRT * v)))
			/ (std::exp(FoRT * v) - 1.0) / (std::exp(FoRT * v) - 1.0)
			- (1.0 * FoRT * Frdy * pK
					* (0.75 * Ko - 0.75 * Ki * std::exp(FoRT * v)))
					/ (std::exp(FoRT * v) - 1.0)
			+ (0.75 * FoRT * FoRT * Frdy * Ki * pK * v * exp(FoRT * v))
					/ (std::exp(FoRT * v) - 1.0);
	dibark *= Q;
	dibark += pK * (v * Frdy * FoRT) * (0.75 * dKi * std::exp(v * FoRT))
			/ (std::exp(v * FoRT) - 1.0);

	double ibarna_j = pNa * (v * Frdy * FoRT)
			* (0.75 * Naj * std::exp(v * FoRT) - 0.75 * Nao)
			/ (std::exp(v * FoRT) - 1.0);
	double dibarna_j = (0.75 * FoRT * FoRT * Frdy * Naj * pNa * v
			* std::exp(FoRT * v)) / (std::exp(FoRT * v) - 1.0)
			- (1.0 * FoRT * Frdy * pNa
					* (0.75 * Nao - 0.75 * Naj * std::exp(FoRT * v)))
					/ (std::exp(FoRT * v) - 1.0)
			+ (FoRT * FoRT * Frdy * pNa * v * std::exp(FoRT * v)
					* (0.75 * Nao - 0.75 * Naj * std::exp(FoRT * v)))
					/ (std::exp(FoRT * v) - 1.0) / (std::exp(FoRT * v) - 1.0);
	dibarna_j *= Q;
	dibarna_j += pNa * (v * Frdy * FoRT) * (0.75 * dNaj * std::exp(v * FoRT))
			/ (std::exp(v * FoRT) - 1.0);

	double ibarna_sl = pNa * (v * Frdy * FoRT)
			* (0.75 * Nasl * std::exp(v * FoRT) - 0.75 * Nao)
			/ (std::exp(v * FoRT) - 1.0);
	double dibarna_sl = (0.75 * FoRT * FoRT * Frdy * Nasl * pNa * v
			* std::exp(FoRT * v)) / (std::exp(FoRT * v) - 1.0)
			- (1.0 * FoRT * Frdy * pNa
					* (0.75 * Nao - 0.75 * Nasl * std::exp(FoRT * v)))
					/ (std::exp(FoRT * v) - 1.0)
			+ (FoRT * FoRT * Frdy * pNa * v * std::exp(FoRT * v)
					* (0.75 * Nao - 0.75 * Nasl * std::exp(FoRT * v)))
					/ (std::exp(FoRT * v) - 1.0) / (std::exp(FoRT * v) - 1.0);
	dibarna_sl *= Q;
	dibarna_sl += pNa * (v * Frdy * FoRT) * (0.75 * dNasl * std::exp(v * FoRT))
			/ (std::exp(v * FoRT) - 1.0);

	aux = std::pow(Q10CaL, Qpow);

	double dI_Ca_junc = (Fjunc_CaL * dibarca_j * d * f * ((1 - fcaBj) + fcaCaj)
			* aux) * 0.45
			+ (Fjunc_CaL * ibarca_j * dd * f * ((1 - fcaBj) + fcaCaj) * aux)
					* 0.45
			+ (Fjunc_CaL * ibarca_j * d * df * ((1 - fcaBj) + fcaCaj) * aux)
					* 0.45
			+ (Fjunc_CaL * ibarca_j * d * f * (-dfcaBj + dfcaCaj) * aux) * 0.45;

	double dI_Ca_sl = (Fsl_CaL * dibarca_sl * d * f * ((1 - fcaBsl) + fcaCaMSL)
			* aux) * 0.45
			+ (Fsl_CaL * ibarca_sl * dd * f * ((1 - fcaBsl) + fcaCaMSL) * aux)
					* 0.45
			+ (Fsl_CaL * ibarca_sl * d * df * ((1 - fcaBsl) + fcaCaMSL) * aux)
					* 0.45
			+ (Fsl_CaL * ibarca_sl * d * f * (-dfcaBsl + dfcaCaMSL) * aux)
					* 0.45;

	double dI_Ca = dI_Ca_junc + dI_Ca_sl;

	double dI_CaK = (dibark * d * f
			* (Fjunc_CaL * (fcaCaj + (1 - fcaBj))
					+ Fsl_CaL * (fcaCaMSL + (1 - fcaBsl))) * aux) * 0.45
			+ (ibark * dd * f
					* (Fjunc_CaL * (fcaCaj + (1 - fcaBj))
							+ Fsl_CaL * (fcaCaMSL + (1 - fcaBsl))) * aux) * 0.45
			+ (ibark * d * df
					* (Fjunc_CaL * (fcaCaj + (1 - fcaBj))
							+ Fsl_CaL * (fcaCaMSL + (1 - fcaBsl))) * aux) * 0.45
			+ (ibark * d * f
					* (Fjunc_CaL * (dfcaCaj + (0 - dfcaBj))
							+ Fsl_CaL * (dfcaCaMSL + (0 - dfcaBsl))) * aux)
					* 0.45;
	double dI_CaNa_junc = (Fjunc_CaL * dibarna_j * d * f
			* ((1 - fcaBj) + fcaCaj) * aux) * 0.45
			+ (Fjunc_CaL * ibarna_j * dd * f * ((1 - fcaBj) + fcaCaj) * aux)
					* 0.45
			+ (Fjunc_CaL * ibarna_j * d * df * ((1 - fcaBj) + fcaCaj) * aux)
					* 0.45
			+ (Fjunc_CaL * ibarna_j * d * f * ((0 - dfcaBj) + dfcaCaj) * aux)
					* 0.45;
	double dI_CaNa_sl = (Fsl_CaL * dibarna_sl * d * f
			* ((1 - fcaBsl) + fcaCaMSL) * aux) * .45
			+ (Fsl_CaL * ibarna_sl * dd * f * ((1 - fcaBsl) + fcaCaMSL) * aux)
					* .45
			+ (Fsl_CaL * ibarna_sl * d * df * ((1 - fcaBsl) + fcaCaMSL) * aux)
					* .45
			+ (Fsl_CaL * ibarna_sl * d * f * ((0 - dfcaBsl) + dfcaCaMSL) * aux)
					* .45;
	double dI_CaNa = dI_CaNa_junc + dI_CaNa_sl;
	double dI_Catot = dI_Ca + dI_CaK + dI_CaNa;

	// I_ncx: Na/Ca Exchanger flux
	double Ka_junc = 1.0 / (1.0 + (Kdact / Caj) * (Kdact / Caj));
	double dKa_junc = -1.0 / (1.0 + (Kdact / Caj) * (Kdact / Caj))
			/ (1.0 + (Kdact / Caj) * (Kdact / Caj))
			* (-2.0 * (Kdact / Caj) * (Kdact / Caj) / Caj * dCaj);
	double Ka_sl = 1.0 / (1.0 + (Kdact / Casl) * (Kdact / Casl));
	double dKa_sl = 1.0 / (1.0 + (Kdact / Casl) * (Kdact / Casl))
			/ (1.0 + (Kdact / Casl) * (Kdact / Casl))
			* (-2.0 * (Kdact / Casl) * (Kdact / Casl) / Casl * dCasl);
	double s1_junc = std::exp(nu * v * FoRT) * Naj * Naj * Naj * Cao;
	double ds1_junc = (nu * FoRT) * s1_junc * Q
			+ std::exp(nu * v * FoRT) * 3 * Naj * Naj * dNaj * Cao;
	double s1_sl = std::exp(nu * v * FoRT) * Nasl * Nasl * Nasl * Cao;
	double ds1_sl = (nu * FoRT) * s1_sl * Q
			+ std::exp(nu * v * FoRT) * 2 * Nasl * Nasl * dNasl * Cao;
	double s2_junc = std::exp((nu - 1) * v * FoRT) * Nao * Nao * Nao * Caj;
	double ds2_junc = ((nu - 1) * FoRT) * s2_junc * Q
			+ std::exp((nu - 1) * v * FoRT) * Nao * Nao * Nao * dCaj;
	double s2_sl = std::exp((nu - 1) * v * FoRT) * Nao * Nao * Nao * Casl;
	double ds2_sl = ((nu - 1) * FoRT) * s2_sl * Q
			+ std::exp((nu - 1) * v * FoRT) * Nao * Nao * Nao * dCasl;

	double s3_junc = KmCai * Nao * Nao * Nao
			* (1 + (Naj / KmNai) * (Naj / KmNai) * (Naj / KmNai))
			+ KmNao * KmNao * KmNao * Caj * (1 + Caj / KmCai)
			+ KmCao * Naj * Naj * Naj + Naj * Naj * Naj * Cao
			+ Nao * Nao * Nao * Caj;
	double ds3_junc = KmCai * Nao * Nao * Nao
			* (1 + 3 * (Naj / KmNai) * (Naj / KmNai) * (dNaj / KmNai))
			+ KmNao * KmNao * KmNao * (dCaj + 2. * Caj * dCaj / KmCai)
			+ KmCao * 3 * Naj * Naj * dNaj + 3 * Naj * Naj * dNaj * Cao
			+ Nao * Nao * Nao * dCaj;
	double s3_sl = KmCai * Nao * Nao * Nao
			* (1 + (Nasl / KmNai) * (Nasl / KmNai) * (Nasl / KmNai))
			+ KmNao * KmNao * KmNao * Casl * (1 + Casl / KmCai)
			+ KmCao * Nasl * Nasl * Nasl + Nasl * Nasl * Nasl * Cao
			+ Nao * Nao * Nao * Casl;
	double ds3_sl = KmCai * Nao * Nao * Nao
			* (1 + 3 * (Nasl / KmNai) * (Nasl / KmNai) * (dNasl / KmNai))
			+ KmNao * KmNao * KmNao * (dCasl + 2 * Casl * dCasl / KmCai)
			+ KmCao * 3 * Nasl * Nasl * dNasl + 3 * Nasl * Nasl * dNasl * Cao
			+ Nao * Nao * Nao * dCasl;

	aux = std::pow(Q10NCX, Qpow);

	double dI_ncx_junc = Fjunc * IbarNCX * aux * dKa_junc * (s1_junc - s2_junc)
			/ s3_junc / (1.0 + ksat * std::exp((nu - 1) * v * FoRT))
			+ Fjunc * IbarNCX * aux * Ka_junc * (ds1_junc - ds2_junc) / s3_junc
					/ (1.0 + ksat * std::exp((nu - 1) * v * FoRT))
			- Fjunc * IbarNCX * aux * Ka_junc * (s1_junc - s2_junc) / s3_junc
					/ s3_junc * ds3_junc
					/ (1.0 + ksat * std::exp((nu - 1) * v * FoRT));
	-Q * Fjunc * IbarNCX * aux * Ka_junc * (s1_junc - s2_junc) / s3_junc
			* (ksat * ((nu - 1) * FoRT) * std::exp((nu - 1) * v * FoRT))
			/ (1.0 + ksat * std::exp((nu - 1) * v * FoRT))
			/ (1.0 + ksat * std::exp((nu - 1) * v * FoRT));

	double dI_ncx_sl = Fsl * IbarNCX * aux * dKa_sl * (s1_sl - s2_sl) / s3_sl
			/ (1 + ksat * std::exp((nu - 1) * v * FoRT))
			+ Fsl * IbarNCX * aux * Ka_sl * (ds1_sl - ds2_sl) / s3_sl
					/ (1 + ksat * std::exp((nu - 1) * v * FoRT))
			- Fsl * IbarNCX * aux * Ka_sl * (s1_sl - s2_sl) / s3_sl / s3_sl
					* ds3_sl / (1 + ksat * std::exp((nu - 1) * v * FoRT))
			- Q * Fsl * IbarNCX * aux * Ka_sl * (s1_sl - s2_sl) / s3_sl
					* (ksat * ((nu - 1) * FoRT) * std::exp((nu - 1) * v * FoRT))
					/ (1 + ksat * std::exp((nu - 1) * v * FoRT))
					/ (1 + ksat * std::exp((nu - 1) * v * FoRT));

	double dI_ncx = dI_ncx_junc + dI_ncx_sl;

	// I_pca: Sarcolemmal Ca Pump Current
	aux = std::pow(Q10SLCaP, Qpow);

	double dI_pca_junc = Fjunc * aux * IbarSLCaP * (1.6 * std::pow(Caj, 0.6))
			/ (std::pow(KmPCa, 1.6) + std::pow(Caj, 1.6))
			- Fjunc * aux * IbarSLCaP * std::pow(Caj, 1.6)
					/ (std::pow(KmPCa, 1.6) + std::pow(Caj, 1.6))
					/ (std::pow(KmPCa, 1.6) + std::pow(Caj, 1.6))
					* (1.6 * std::pow(Caj, 0.6));
	dI_pca_junc *= dCaj;

	double dI_pca_sl = Fsl * aux * IbarSLCaP * (1.6 * std::pow(Casl, 0.6))
			/ (std::pow(KmPCa, 1.6) + std::pow(Casl, 1.6))
			- Fsl * aux * IbarSLCaP * std::pow(Casl, 1.6)
					/ (std::pow(KmPCa, 1.6) + std::pow(Casl, 1.6))
					/ (std::pow(KmPCa, 1.6) + std::pow(Casl, 1.6))
					* (1.6 * std::pow(Casl, 0.6));
	dI_pca_sl *= dCasl;
	double dI_pca = dI_pca_junc + dI_pca_sl;
	// I_cabk: Ca Background Current
	double dI_cabk_junc = Fjunc * GCaB * Q;
	double dI_cabk_sl = Fsl * GCaB * Q;
	double dI_cabk = dI_cabk_junc + dI_cabk_sl;

	// Sodium Concentrations
	double dI_Na_tot_junc = dI_Na_junc + dI_nabk_junc + 3 * dI_ncx_junc
			+ 3 * dI_nak_junc + dI_CaNa_junc + dI_NaL_junc;   // [uA/uF]
	double dI_Na_tot_sl = dI_Na_sl + dI_nabk_sl + 3 * dI_ncx_sl + 3 * dI_nak_sl
			+ dI_CaNa_sl + dI_NaL_sl;   // [uA/uF]
	double dI_Na_tot_sl2 = 3 * dI_ncx_sl + 3 * dI_nak_sl + dI_CaNa_sl; // [uA/uF]
	double dI_Na_tot_junc2 = 3 * dI_ncx_junc + 3 * dI_nak_junc + dI_CaNa_junc; // [uA/uF]

	// Potassium Concentration
	double dI_K_tot = dI_to + dI_kr + dI_ks + dI_ki - 2 * dI_nak + dI_CaK
			+ dI_kp + dI_kur;     // [uA/uF] //SVP: added IKur

	// Calcium Concentrations
	double dI_Ca_tot_junc = dI_Ca_junc + dI_cabk_junc + dI_pca_junc
			- 2 * dI_ncx_junc;                   // [uA/uF]
	double dI_Ca_tot_sl = dI_Ca_sl + dI_cabk_sl + dI_pca_sl - 2 * dI_ncx_sl; // [uA/uF]

	double dI_Na_tot = dI_Na_tot_junc + dI_Na_tot_sl;          // [uA/uF]
	double dI_Cl_tot = dI_ClCa + dI_Clbk + dI_ClCFTR;                 // [uA/uF]
	double dI_Ca_tot = dI_Ca_tot_junc + dI_Ca_tot_sl;
	double dI_tot = dI_Na_tot + dI_Cl_tot + dI_Ca_tot + dI_K_tot;

	return dI_tot;
}

void Kharche11::updateVariables(std::vector<double>& variables,
		double appliedCurrent, double dt)
{
	nernst_potentials(variables);

	// Membrane Currents
	double& v = variables[0];
	double& m = variables[1];
	double& h = variables[2];
	double& j = variables[3];
	double& d = variables[4];
	double& f = variables[5];
	double& fcaBj = variables[6];
	double& fcaBsl = variables[7];
	double& xtof = variables[8];
	double& ytof = variables[9];
	double& xkr = variables[10];
	double& xks = variables[11];
	double& RyRr = variables[12];
	double& RyRo = variables[13];
	double& RyRi = variables[14];
	double& NaBj = variables[15];
	double& NaBsl = variables[16];
	double& TnCL = variables[17];
	double& TnCHc = variables[18];
	double& TnCHm = variables[19];
	double& CaM = variables[20];
	double& Myoc = variables[21];
	double& Myom = variables[22];
	double& SRB = variables[23];
	double& SLLj = variables[24];
	double& SLLsl = variables[25];
	double& SLHj = variables[26];
	double& SLHsl = variables[27];
	double& Csqnb = variables[28];
	double& Ca_sr = variables[29];
	double& Naj = variables[30];
	double& Nasl = variables[31];
	double& Nai = variables[32];
	double& Ki = variables[33];
	double& Caj = variables[34];
	double& Casl = variables[35];
	double& Cai = variables[36];
	double& C1 = variables[37];
	double& C2 = variables[38];
	double& C3 = variables[39];
	double& C4 = variables[40];
	double& C5 = variables[41];
	double& C6 = variables[42];
	double& C7 = variables[43];
	double& C8 = variables[44];
	double& C9 = variables[45];
	double& C10 = variables[46];
	double& C11 = variables[47];
	double& C12 = variables[48];
	double& C13 = variables[49];
	double& C14 = variables[50];
	double& C15 = variables[51];
	double& O1 = variables[52];
	double& rkur = variables[53];
	double& skur = variables[54];
	double& Inal1 = variables[55];
	double& Inal2 = variables[56];

	double aux = (1. + std::exp(-(56.86 + v) / 9.03));
	double mss = 1. / (aux * aux);

	aux = (v + 45.79) / 15.54;
	double taum = 0.1292 * std::exp(-aux * aux);
	aux = (v - 4.823) / 51.12;
	taum += 0.06487 * std::exp(-aux * aux);

	double ah = (v < -40.0) ? 0.057 * std::exp(-(v + 80) / 6.8) : 0.0;
	double bh =
			(v < -40.0) ?
					2.7 * std::exp(0.079 * v)
							+ 3.1 * 1e5 * std::exp(0.3485 * v) :
					0.77 / (0.13 * (1.0 + std::exp(-(v + 10.66) / 11.1)));

	double tauh = 1.0 / (ah + bh);

	aux = 1.0 + std::exp((v + 71.55) / 7.43);
	double hss = 1.0 / (aux * aux);

	//    double  aj = (v >= -40) * (0)
	//                +(v < -40) * (((-2.5428 * 10^4*std::exp(0.2444*v) - 6.948*10^-6 * std::exp(-0.04391*v)) * (v + 37.78)) /
	//                         (1 + std::exp( 0.311 * (v + 79.23) )));
	double aj =
			(v < -40.0) ?
					(-2.5428 * 1e4 * std::exp(0.2444 * v)
							- 6.948 * 1e-6 * std::exp(-0.04391 * v))
							* (v + 37.78)
							/ (1.0 + std::exp(0.311 * (v + 79.23))) :
					0.0;
	//    double bj = (v >= -40) * ((0.6 * std::exp( 0.057 * v)) / (1 + std::exp( -0.1 * (v + 32) )))
	//       + (v < -40) * ((0.02424 * std::exp( -0.01052 * v )) / (1 + std::exp( -0.1378 * (v + 40.14) )));
	double bj =
			(v < -40.0) ?
					((0.02424 * std::exp(-0.01052 * v))
							/ (1 + std::exp(-0.1378 * (v + 40.14)))) :
					((0.6 * std::exp(0.057 * v))
							/ (1 + std::exp(-0.1 * (v + 32))));
	double tauj = 1.0 / (aj + bj);
	aux = 1.0 + std::exp((v + 71.55) / 7.43);
	double jss = 1.0 / (aux * aux);

	// TODO:: USE RL
	// sm = M_INF-(M_INF-sm)*std::exp(-dt/TAU_M);

	//ydot(1) = (mss - y(1)) / taum;
	//m += dt * (mss - m) / taum;
	 m = mss - (mss - m) * std::exp(-dt / taum);
	//ydot(2) = (hss - y(2)) / tauh;
	h = hss - (hss - h) * std::exp(-dt / tauh);
	//h += dt * (hss - h) / tauh;
	//ydot(3) = (jss - y(3)) / tauj;
	j = jss - (jss - j) * std::exp(-dt / tauj);
	//j += dt * (jss - j) / tauj;

	I_Na_junc = Fjunc * GNa * m * m * m * h * j * (v - ena_junc);
	I_Na_sl = Fsl * GNa * m * m * m * h * j * (v - ena_sl);
	I_Na = I_Na_junc + I_Na_sl;

	// Late I_Na
	double aml = 0.32 * (v + 47.13) / (1 - std::exp(-0.1 * (v + 47.13)));
	double bml = 0.08 * std::exp(-v / 11);
	double hlinf = 1 / (1 + std::exp((v + 91) / 6.1));
	double tauhl = 600.0;
	//TODO use RL
	//ydot(60) = aml*(1-y(60))-bml*y(60);
	double tau_Inal1 = 1.0 / (aml + bml);
	double Inal1_inf = aml * tau_Inal1;
	Inal1 = Inal1_inf - (Inal1_inf - Inal1) * std::exp(-dt / tau_Inal1);
	//Inal1 +=  dt * ( aml*(1-Inal1)-bml*Inal1 );

	//ydot(61) = (hlinf-y(61))/tauhl;
	Inal2 = hlinf - (hlinf - Inal2) * std::exp(-dt / tauhl);
	//Inal2 +=  dt * ( hlinf - Inal2 ) / tauhl;

	I_NaL_junc = Fjunc * GNaL * Inal1 * Inal1 * Inal1 * Inal2 * (v - ena_junc);
	I_NaL_sl = Fsl * GNaL * Inal1 * Inal1 * Inal1 * Inal2 * (v - ena_sl);
	I_NaL = I_NaL_junc + I_NaL_sl;
	// No idea about this
	//if t<9050
	//    ydot(62)=0;
	//else
	//    ydot(62)=I_NaL;
	//end
	// I_nabk: Na Background Current
	I_nabk_junc = Fjunc * GNaB * (v - ena_junc);
	I_nabk_sl = Fsl * GNaB * (v - ena_sl);
	I_nabk = I_nabk_junc + I_nabk_sl;

	// I_nak: Na/K Pump Current
	double sigma = (std::exp(Nao / 67.3) - 1) / 7.0;
	double fnak = 1.0
			/ (1.0 + 0.1245 * std::exp(-0.1 * v * FoRT)
					+ 0.0365 * sigma * std::exp(-v * FoRT));
	aux = KmNaip / Naj;
	I_nak_junc = Fjunc * IbarNaK * fnak * Ko / (1.0 + aux * aux * aux * aux)
			/ (Ko + KmKo);
	aux = KmNaip / Nasl;
	I_nak_sl = Fsl * IbarNaK * fnak * Ko / (1.0 + aux * aux * aux * aux)
			/ (Ko + KmKo);
	I_nak = I_nak_junc + I_nak_sl;

	// I_kr: Rapidly Activating K Current
	double gkr = 0.035 * std::sqrt(Ko / 5.4);
	double xrss = 1.0 / (1.0 + std::exp(-(v + 10.0) / 5.0));
	double tauxr = 550.0 / (1.0 + std::exp((-22.0 - v) / 9.0)) * 6.0
			/ (1.0 + std::exp((v - (-11.0)) / 9.0))
			+ 230.0 / (1.0 + std::exp((v - (-40.0)) / 20.0));
	// TODO: Use RL
	//ydot(12) = (xrss-y(12))/tauxr;
	xkr = xrss - (xrss - xkr) * std::exp(-dt / tauxr);
	//xkr += dt * (xrss - xkr) / tauxr ;
	double rkr = 1.0 / (1.0 + std::exp((v + 74.0) / 24.0));
	I_kr = gkr * xkr * rkr * (v - ek);

	// I_ks: Slowly Activating K Current

	double eks = (1 / FoRT) * std::log((Ko + pNaK * Nao) / (Ki + pNaK * Nai));

	if (!markov_iks)
	{
		double gks_junc = 1.0 * (1.0 + 1.0 * AF + 2.0 * ISO) * 0.0035 * 1.0;
		double gks_sl = 1.0 * (1.0 + 1.0 * AF + 2.0 * ISO) * 0.0035 * 1.0; // Fra
		double xsss = 1.0 / (1.0 + std::exp(-(v + 40.0 * ISO + 3.8) / 14.25)); // fitting Fra
		double tauxs = 990.1
				/ (1 + std::exp(-(v + 40.0 * ISO + 2.436) / 14.12));
		// TODO: Use RL
		//ydot(13) = (xsss-y(13))/tauxs;
		xks = xsss - (xsss - xks) * std::exp(-dt / tauxs);
		//xks += dt * (xsss - xks) / tauxs ;

		I_ks_junc = Fjunc * gks_junc * xks * xks * (v - eks);
		I_ks_sl = Fsl * gks_sl * xks * xks * (v - eks);
		I_ks = I_ks_junc + I_ks_sl;
	}
	else
	{
		double gks_junc = 1 * 0.0065;
		double gks_sl = 1 * 0.0065; //FRA
		double alpha = 3.98e-4 * std::exp(3.61e-1 * v * FoRT);
		double beta = 5.74e-5 * std::exp(-9.23e-2 * v * FoRT);
		double gamma = 3.41e-3 * std::exp(8.68e-1 * v * FoRT);
		double delta = 1.2e-3 * std::exp(-3.3e-1 * v * FoRT);
		double teta = 6.47e-3;
		double eta = 1.25e-2 * std::exp(-4.81e-1 * v * FoRT);
		double psi = 6.33e-3 * std::exp(1.27 * v * FoRT);
		double omega = 4.91e-3 * std::exp(-6.79e-1 * v * FoRT);

		//ydot(42)=-4*alpha*C1+beta*C2;
		//ydot(43)=4*alpha*C1-(beta+gamma+3*alpha)*C2+2*beta*C3;
		//ydot(44)=3*alpha*C2-(2*beta+2*gamma+2*alpha)*C3+3*beta*C4;
		//ydot(45)=2*alpha*C3-(3*beta+3*gamma+alpha)*C4+4*beta*C5;
		//ydot(46)=1*alpha*C3-(4*beta+4*gamma)*C5+delta*C9;
		//ydot(47)=gamma*C2-(delta+3*alpha)*C6+beta*C7;
		//ydot(48)=2*gamma*C3+3*alpha*C6-(delta+beta+2*alpha+gamma)*C7+2*beta*C8+2*delta*C10;
		//ydot(49)=3*gamma*C4+2*alpha*C7-(delta+2*beta+1*alpha+2*gamma)*C8+3*beta*C9+2*delta*C11;
		//ydot(50)=4*gamma*C5+1*alpha*C8-(delta+3*beta+0*alpha+3*gamma)*C9+2*delta*C12;
		//ydot(51)=1*gamma*C7-(2*delta+2*alpha)*C10+beta*C11;
		//ydot(52)=2*gamma*C8+2*alpha*C10-(2*delta+beta+1*alpha+gamma)*C11+2*beta*C12+3*delta*C13;
		//ydot(53)=3*gamma*C9+1*alpha*C11-(2*delta+2*beta+2*gamma)*C12+3*delta*C14;
		//ydot(54)=1*gamma*C11-(3*delta+1*alpha)*C13+beta*C14;
		//ydot(55)=2*gamma*C12+1*alpha*C13-(3*delta+1*beta+1*gamma)*C14+4*delta*C15;
		//ydot(56)=1*gamma*C14-(4*delta+teta)*C15+eta*O1;
		//O2=1-(C1+C2+C3+C4+C5+C6+C8+C7+C9+C10+C11+C12+C13+C14+C15+O1);
		//ydot(57)=1*teta*C15-(eta+psi)*O1+omega*O2;

		//dC1 =-4*alpha*C1+beta*C2;
		//dC2 =4*alpha*C1-(beta+gamma+3*alpha)*C2+2*beta*C3;
		//dC3 =3*alpha*C2-(2*beta+2*gamma+2*alpha)*C3+3*beta*C4;
		//dC4 =2*alpha*C3-(3*beta+3*gamma+alpha)*C4+4*beta*C5;
		//dC5 =1*alpha*C3-(4*beta+4*gamma)*C5+delta*C9;
		//dC6 =gamma*C2-(delta+3*alpha)*C6+beta*C7;
		//dC7 =2*gamma*C3+3*alpha*C6-(delta+beta+2*alpha+gamma)*C7+2*beta*C8+2*delta*C10;
		//dC8 =3*gamma*C4+2*alpha*C7-(delta+2*beta+1*alpha+2*gamma)*C8+3*beta*C9+2*delta*C11;
		//dC9 =4*gamma*C5+1*alpha*C8-(delta+3*beta+0*alpha+3*gamma)*C9+2*delta*C12;
		//dC10=1*gamma*C7-(2*delta+2*alpha)*C10+beta*C11;
		//dC11=2*gamma*C8+2*alpha*C10-(2*delta+beta+1*alpha+gamma)*C11+2*beta*C12+3*delta*C13;
		//dC12=3*gamma*C9+1*alpha*C11-(2*delta+2*beta+2*gamma)*C12+3*delta*C14;
		//dC13=1*gamma*C11-(3*delta+1*alpha)*C13+beta*C14;
		//dC14=2*gamma*C12+1*alpha*C13-(3*delta+1*beta+1*gamma)*C14+4*delta*C15;
		//dC15=1*gamma*C14-(4*delta+teta)*C15+eta*O1;
		//double O2=1-(C1+C2+C3+C4+C5+C6+C8+C7+C9+C10+C11+C12+C13+C14+C15+O1);
		//dO1 =1*teta*C15-(eta+psi)*O1+omega*O2;

		C1 += dt * (-4 * alpha * C1 + beta * C2);
		C2 += dt
				* (4 * alpha * C1 - (beta + gamma + 3 * alpha) * C2
						+ 2 * beta * C3);
		C3 += dt
				* (3 * alpha * C2 - (2 * beta + 2 * gamma + 2 * alpha) * C3
						+ 3 * beta * C4);
		C4 += dt
				* (2 * alpha * C3 - (3 * beta + 3 * gamma + alpha) * C4
						+ 4 * beta * C5);
		C5 += dt * (1 * alpha * C3 - (4 * beta + 4 * gamma) * C5 + delta * C9);
		C6 += dt * (gamma * C2 - (delta + 3 * alpha) * C6 + beta * C7);
		C7 += dt
				* (2 * gamma * C3 + 3 * alpha * C6
						- (delta + beta + 2 * alpha + gamma) * C7
						+ 2 * beta * C8 + 2 * delta * C10);
		C8 += dt
				* (3 * gamma * C4 + 2 * alpha * C7
						- (delta + 2 * beta + 1 * alpha + 2 * gamma) * C8
						+ 3 * beta * C9 + 2 * delta * C11);
		C9 += dt
				* (4 * gamma * C5 + 1 * alpha * C8
						- (delta + 3 * beta + 0 * alpha + 3 * gamma) * C9
						+ 2 * delta * C12);
		C10 += dt
				* (1 * gamma * C7 - (2 * delta + 2 * alpha) * C10 + beta * C11);
		C11 += dt
				* (2 * gamma * C8 + 2 * alpha * C10
						- (2 * delta + beta + 1 * alpha + gamma) * C11
						+ 2 * beta * C12 + 3 * delta * C13);
		C12 += dt
				* (3 * gamma * C9 + 1 * alpha * C11
						- (2 * delta + 2 * beta + 2 * gamma) * C12
						+ 3 * delta * C14);
		C13 +=
				dt
						* (1 * gamma * C11 - (3 * delta + 1 * alpha) * C13
								+ beta * C14);
		C14 += dt
				* (2 * gamma * C12 + 1 * alpha * C13
						- (3 * delta + 1 * beta + 1 * gamma) * C14
						+ 4 * delta * C15);
		C15 += dt * (1 * gamma * C14 - (4 * delta + teta) * C15 + eta * O1);
		double O2 = 1
				- (C1 + C2 + C3 + C4 + C5 + C6 + C8 + C7 + C9 + C10 + C11 + C12
						+ C13 + C14 + C15 + O1);
		O1 += dt * (1 * teta * C15 - (eta + psi) * O1 + omega * O2);
		I_ks_junc = Fjunc * gks_junc * (O1 + O2) * (v - eks);
		I_ks_sl = Fsl * gks_sl * (O1 + O2) * (v - eks);
		I_ks = I_ks_junc + I_ks_sl;
	}

	// I_kp: Plateau K current
	double kp_kp = 1.0 / (1.0 + std::exp(7.488 - v / 5.98));
	I_kp_junc = Fjunc * gkp * kp_kp * (v - ek);
	I_kp_sl = Fsl * gkp * kp_kp * (v - ek);
	I_kp = I_kp_junc + I_kp_sl;

	// I_to: Transient Outward K Current (slow and fast components)
	// modified for human myocytes

	// 11/12/09; changed Itof to that from maleckar/giles/2009; removed I_tos
	// atrium
	// equations for activation;
	double xtoss = (1.0 / (1 + std::exp(-(v + 1.0) / 11.0)));
	double tauxtof = 3.5 * std::exp(-((v / 30.0) * (v / 30.0))) + 1.5;

	// equations for inactivation;
	double ytoss = (1.0 / (1 + std::exp((v + 40.5) / 11.5)));
	double tauytof = 25.635
			* std::exp(-(((v + 52.45) / 15.8827) * ((v + 52.45) / 15.8827)))
			+ 24.14; //14.14

	// TODO: Use RL
	//ydot(10) = (xtoss-y(10))/tauxtof;
	xtof = xtoss - (xtoss - xtof) * std::exp(-dt / tauxtof);
	//xtof += dt * ( xtoss - xtof ) / tauxtof ;
	//ydot(11) = (ytoss-y(11))/tauytof;
	ytof = ytoss - (ytoss - ytof) * std::exp(-dt / tauytof);
	//ytof += dt *( ytoss - ytof )  / tauytof;
	I_tof = 1.0 * GtoFast * xtof * ytof * (v - ek);
	I_to = 1.0 * I_tof;

	// I_kur: Ultra rapid delayed rectifier Outward K Current
	// Equation for IKur; from Maleckar et al. 2009 - EG
	// atrium
	// equations for activation;
	double RA = 0; // Right Atrium
	double Gkur = 1 * (1.0 - 0.5 * AF) * (1 + 2 * ISO) * 0.045 * (1 + 0.2 * RA); // nS/pF maleckar 0.045
	double xkurss = (1.0 / (1 + std::exp((v + 6) / -8.6)));
	double tauxkur = 9 / (1 + std::exp((v + 5) / 12.0)) + 0.5;

	// equations for inactivation;
	double ykurss = (1.0 / (1 + std::exp((v + 7.5) / 10)));
	double tauykur = 590 / (1 + std::exp((v + 60) / 10.0)) + 3050;

	// TODO: Use RL
	//ydot(58) = (xkurss-y(58))/tauxkur;
	 rkur = xkurss - (xkurss - rkur) * std::exp(-dt / tauxkur);
	//rkur += dt * ( xkurss - rkur )  / tauxkur;
	//ydot(59) = (ykurss-y(59))/tauykur;
	  skur = ykurss - (ykurss - skur) * std::exp(-dt / tauykur);
	//skur += dt * ( ykurss - skur )  / tauykur;

	I_kur = 1 * Gkur * rkur * skur * (v - ek);

	// I_ki: Time-Independent K Current
	double aki = 1.02 / (1 + std::exp(0.2385 * (v - ek - 59.215)));
	double bki = (0.49124 * std::exp(0.08032 * (v + 5.476 - ek))
			+ std::exp(0.06175 * (v - ek - 594.31)))
			/ (1.0 + std::exp(-0.5143 * (v - ek + 4.753)));
	double kiss = aki / (aki + bki);

	// I_ki =1* 0.35*sqrt(Ko/5.4)*kiss*(v-ek);
	// SVP 11/11/09
	// multiplieD IK1 by 0.15 to scale it to single cell isolated atrial cell
	// resting potential
	I_ki = (1.0 + 1.0 * AF) * 0.0525 * std::sqrt(Ko / 5.4) * kiss * (v - ek);

	// I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
	I_ClCa_junc = Fjunc * GClCa / (1 + KdClCa / Caj) * (v - ecl);
	I_ClCa_sl = Fsl * GClCa / (1 + KdClCa / Casl) * (v - ecl);
	I_ClCa = I_ClCa_junc + I_ClCa_sl;
	I_Clbk = GClB * (v - ecl);

	double GClCFTR = 0; // 4.9e-3*ISO;     // [mS/uF]
	I_ClCFTR = GClCFTR * (v - ecl);

	// I_Ca: L-type Calcium Current
	double dss = 1.0 / (1.0 + std::exp(-(v + 3 * ISO + 9) / 6)); //in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
	double taud = 1.0 * dss * (1 - std::exp(-(v + 3 * ISO + 9) / 6))
			/ (0.035 * (v + 3 * ISO + 9));
	double fss = 1.0 / (1.0 + std::exp((v + 3 * ISO + 30) / 7))
			+ 0.2 / (1 + std::exp((50 - v - 3 * ISO) / 20)); // in Maleckar v1/2=-27.4 S=7.1 (mV); Courtemanche v1/2=-28 S=6.9 (mV)
	double tauf = 1.0
			/ (0.0197
					* std::exp(
							-(0.0337 * (v + 3 * ISO + 25))
									* (0.0337 * (v + 3 * ISO + 25))) + 0.02);

	// TODO: Use RL
	//ydot(4) = (dss-d)/taud;
	d = dss - (dss - d) * std::exp(-dt / taud);
	//d += dt * ( dss - d )   / taud;
	//ydot(5) = (fss-f)/tauf;
	f = fss - (fss - f) * std::exp(-dt / tauf);
	//f += dt *  ( fss - f )  / tauf;
	//ydot(6) = 1.7*Caj*(1-fcaBj)-1*11.9e-3*fcaBj; // fCa_junc   koff!!!!!!!!
	fcaBj += dt * (1.7 * Caj * (1 - fcaBj) - 1 * 11.9e-3 * fcaBj); // fCa_junc   koff!!!!!!!!
	//ydot(7) = 1.7*Casl*(1-fcaBsl)-1*11.9e-3*fcaBsl; // fCa_sl
	fcaBsl += dt * (1.7 * Casl * (1 - fcaBsl) - 1 * 11.9e-3 * fcaBsl);

	// ////////////////////////////////////////////////////////
	// ////////////////////////////////////////////////////////
	double fcaCaMSL = 0.1 / (1.0 + (0.01 / Casl));
	double fcaCaj = 0.1 / (1.0 + (0.01 / Caj));
	fcaCaMSL = 0.0;
	fcaCaj = 0.0;
	double ibarca_j = pCa * 4 * (v * Frdy * FoRT)
			* (0.341 * Caj * std::exp(2 * v * FoRT) - 0.341 * Cao)
			/ (std::exp(2 * v * FoRT) - 1.0);
	double ibarca_sl = pCa * 4 * (v * Frdy * FoRT)
			* (0.341 * Casl * std::exp(2 * v * FoRT) - 0.341 * Cao)
			/ (std::exp(2 * v * FoRT) - 1.0);
	double ibark = pK * (v * Frdy * FoRT)
			* (0.75 * Ki * std::exp(v * FoRT) - 0.75 * Ko)
			/ (std::exp(v * FoRT) - 1.0);
	double ibarna_j = pNa * (v * Frdy * FoRT)
			* (0.75 * Naj * std::exp(v * FoRT) - 0.75 * Nao)
			/ (std::exp(v * FoRT) - 1.0);
	double ibarna_sl = pNa * (v * Frdy * FoRT)
			* (0.75 * Nasl * std::exp(v * FoRT) - 0.75 * Nao)
			/ (std::exp(v * FoRT) - 1.0);
	aux = std::pow(Q10CaL, Qpow);
	I_Ca_junc = (Fjunc_CaL * ibarca_j * d * f * ((1 - fcaBj) + fcaCaj) * aux)
			* 0.45;
	I_Ca_sl = (Fsl_CaL * ibarca_sl * d * f * ((1 - fcaBsl) + fcaCaMSL) * aux)
			* 0.45;
	I_Ca = I_Ca_junc + I_Ca_sl;
	I_CaK = (ibark * d * f
			* (Fjunc_CaL * (fcaCaj + (1 - fcaBj))
					+ Fsl_CaL * (fcaCaMSL + (1 - fcaBsl))) * aux) * 0.45;
	I_CaNa_junc = (Fjunc_CaL * ibarna_j * d * f * ((1 - fcaBj) + fcaCaj) * aux)
			* 0.45;
	I_CaNa_sl = (Fsl_CaL * ibarna_sl * d * f * ((1 - fcaBsl) + fcaCaMSL) * aux)
			* .45;
	I_CaNa = I_CaNa_junc + I_CaNa_sl;
	I_Catot = I_Ca + I_CaK + I_CaNa;

	// I_ncx: Na/Ca Exchanger flux
	double Ka_junc = 1.0 / (1.0 + (Kdact / Caj) * (Kdact / Caj));
	double Ka_sl = 1.0 / (1.0 + (Kdact / Casl) * (Kdact / Casl));
	double s1_junc = std::exp(nu * v * FoRT) * Naj * Naj * Naj * Cao;
	double s1_sl = std::exp(nu * v * FoRT) * Nasl * Nasl * Nasl * Cao;
	double s2_junc = std::exp((nu - 1) * v * FoRT) * Nao * Nao * Nao * Caj;
	double s3_junc = KmCai * Nao * Nao * Nao
			* (1 + (Naj / KmNai) * (Naj / KmNai) * (Naj / KmNai))
			+ KmNao * KmNao * KmNao * Caj * (1 + Caj / KmCai)
			+ KmCao * Naj * Naj * Naj + Naj * Naj * Naj * Cao
			+ Nao * Nao * Nao * Caj;
	double s2_sl = std::exp((nu - 1) * v * FoRT) * Nao * Nao * Nao * Casl;
	double s3_sl = KmCai * Nao * Nao * Nao
			* (1 + (Nasl / KmNai) * (Nasl / KmNai) * (Nasl / KmNai))
			+ KmNao * KmNao * KmNao * Casl * (1 + Casl / KmCai)
			+ KmCao * Nasl * Nasl * Nasl + Nasl * Nasl * Nasl * Cao
			+ Nao * Nao * Nao * Casl;

	aux = std::pow(Q10NCX, Qpow);
	I_ncx_junc = Fjunc * IbarNCX * aux * Ka_junc * (s1_junc - s2_junc) / s3_junc
			/ (1.0 + ksat * std::exp((nu - 1) * v * FoRT));
	I_ncx_sl = Fsl * IbarNCX * aux * Ka_sl * (s1_sl - s2_sl) / s3_sl
			/ (1 + ksat * std::exp((nu - 1) * v * FoRT));
	I_ncx = I_ncx_junc + I_ncx_sl;

	// I_pca: Sarcolemmal Ca Pump Current
	aux = std::pow(Q10SLCaP, Qpow);
	I_pca_junc = Fjunc * aux * IbarSLCaP * std::pow(Caj, 1.6)
			/ (std::pow(KmPCa, 1.6) + std::pow(Caj, 1.6));
	I_pca_sl = Fsl * aux * IbarSLCaP * std::pow(Casl, 1.6)
			/ (std::pow(KmPCa, 1.6) + std::pow(Casl, 1.6));
	I_pca = I_pca_junc + I_pca_sl;

	// I_cabk: Ca Background Current
	I_cabk_junc = Fjunc * GCaB * (v - eca_junc);
	I_cabk_sl = Fsl * GCaB * (v - eca_sl);
	I_cabk = I_cabk_junc + I_cabk_sl;

	// SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
	double MaxSR = 15;
	double MinSR = 1;
	double kCaSR = MaxSR
			- (MaxSR - MinSR) / (1 + std::pow(ec50SR / Ca_sr, 2.5));
	double koSRCa = koCa / kCaSR;
	double kiSRCa = kiCa * kCaSR;
	double RI = 1 - RyRr - RyRo - RyRi;

	//ydot(14) = (kim*RI-kiSRCa*Caj*RyRr)-(koSRCa*Caj*Caj*RyRr-kom*RyRo);   // R
	//ydot(15) = (koSRCa*Caj*Caj*RyRr-kom*RyRo)-(kiSRCa*Caj*RyRo-kim*RyRi);// O
	//ydot(16) = (kiSRCa*Caj*RyRo-kim*RyRi)-(kom*RyRi-koSRCa*Caj*Caj*RI);   // I
	RyRr += dt
			* ((kim * RI - kiSRCa * Caj * RyRr)
					- (koSRCa * Caj * Caj * RyRr - kom * RyRo));   // R
	RyRo += dt
			* ((koSRCa * Caj * Caj * RyRr - kom * RyRo)
					- (kiSRCa * Caj * RyRo - kim * RyRi));   // O
	RyRi += dt
			* ((kiSRCa * Caj * RyRo - kim * RyRi)
					- (kom * RyRi - koSRCa * Caj * Caj * RI));   // I

	J_SRCarel = ks * RyRo * (Ca_sr - Caj);          // [mM/ms]

	aux = std::pow(Q10SRCaP, Qpow);
	J_serca = 1.0 * aux * Vmax_SRCaP
			* (std::pow((Cai / Kmf), hillSRCaP)
					- std::pow((Ca_sr / Kmr), hillSRCaP))
			/ (1.0 + std::pow((Cai / Kmf), hillSRCaP)
					+ std::pow((Ca_sr / Kmr), hillSRCaP));
	J_SRleak = (1.0 + 0.25 * AF) * 5.348e-6 * (Ca_sr - Caj);        //   [mM/ms]

	// Sodium and Calcium Buffering
	//ydot(17) = kon_na*Naj*(Bmax_Naj-NaBj)-koff_na*NaBj;        // NaBj      [mM/ms]
	//ydot(18) = kon_na*Nasl*(Bmax_Nasl-NaBsl)-koff_na*NaBsl;       // NaBsl     [mM/ms]
	double dNaBj = (kon_na * Naj * (Bmax_Naj - NaBj) - koff_na * NaBj);
	NaBj += dt * dNaBj;        // NaBj      [mM/ms]
	double dNaBsl = (kon_na * Nasl * (Bmax_Nasl - NaBsl) - koff_na * NaBsl);
	NaBsl += dt * dNaBsl;       // NaBsl     [mM/ms]

	// Cytosolic Ca Buffers
	double dTnCL = kon_tncl * Cai * (Bmax_TnClow - TnCL) - koff_tncl * TnCL; // TnCL      [mM/ms]
	TnCL += dt * dTnCL;
	double dTnCHc = kon_tnchca * Cai * (Bmax_TnChigh - TnCHc - TnCHm)
			- koff_tnchca * TnCHc; // TnCHc     [mM/ms]
	TnCHc += dt * dTnCHc;
	double dTnCHm = kon_tnchmg * Mgi * (Bmax_TnChigh - TnCHc - TnCHm)
			- koff_tnchmg * TnCHm;   // TnCHm     [mM/ms]
	TnCHm += dt * dTnCHm;
	double dCaM = kon_cam * Cai * (Bmax_CaM - CaM) - koff_cam * CaM; // CaM       [mM/ms]
	CaM += dt * dCaM;
	double dMyoc = kon_myoca * Cai * (Bmax_myosin - Myoc - Myom)
			- koff_myoca * Myoc;    // Myosin_ca [mM/ms]
	Myoc += dt * dMyoc;
	double dMyom = kon_myomg * Mgi * (Bmax_myosin - Myoc - Myom)
			- koff_myomg * Myom;      // Myosin_mg [mM/ms]
	Myom += dt * dMyom;
	double dSRB = kon_sr * Cai * (Bmax_SR - SRB) - koff_sr * SRB; // SRB       [mM/ms]
	SRB += dt * dSRB;
	J_CaB_cytosol = dTnCL + dTnCHc + dTnCHm + dCaM + dMyoc + dMyom + dSRB;

	// Junctional and SL Ca Buffers
	double dSLLj = kon_sll * Caj * (Bmax_SLlowj - SLLj) - koff_sll * SLLj; // SLLj      [mM/ms]
	SLLj += dt * dSLLj;
	double dSLLsl = kon_sll * Casl * (Bmax_SLlowsl - SLLsl) - koff_sll * SLLsl; // SLLsl     [mM/ms]
	SLLsl += dt * dSLLsl;
	double dSLHj = kon_slh * Caj * (Bmax_SLhighj - SLHj) - koff_slh * SLHj; // SLHj      [mM/ms]
	SLHj += dt * dSLHj;
	double dSLHsl = kon_slh * Casl * (Bmax_SLhighsl - SLHsl) - koff_slh * SLHsl; // SLHsl     [mM/ms]
	SLHsl += dt * dSLHsl;
	J_CaB_junction = dSLLj + dSLHj;
	J_CaB_sl = dSLLsl + dSLHsl;

	//// Ion concentrations
	// SR Ca Concentrations
	double dCsqnb = kon_csqn * Ca_sr * (Bmax_Csqn - Csqnb) - koff_csqn * Csqnb; // Csqn      [mM/ms]
	Csqnb += dt * dCsqnb;
	//Csqnb = (Csqnb + dt * kon_csqn * Ca_sr * Bmax_Csqn)
	//		/ (1.0 + dt * (kon_csqn * Ca_sr + koff_csqn));
	double dCa_sr = J_serca - (J_SRleak * Vmyo / Vsr + J_SRCarel) - dCsqnb; // Ca_sr     [mM/ms] //Ratio 3 leak current
	Ca_sr += dt * dCa_sr;

//    J_SRCarel = ks*RyRo*(Ca_sr-Caj);          // [mM/ms]
//
//    aux = std::pow(Q10SRCaP,Qpow);
//    J_serca = 1.0*aux*Vmax_SRCaP*(std::pow( (Cai/Kmf), hillSRCaP)-std::pow( (Ca_sr/Kmr) , hillSRCaP) )
//        /(1.0 + std::pow( (Cai/Kmf), hillSRCaP) + std::pow( (Ca_sr/Kmr), hillSRCaP) );
//    J_SRleak = (1.0+0.25*AF)*5.348e-6*(Ca_sr-Caj);           //   [mM/ms]

//	Ca_sr = (Ca_sr
//			+ dt
//					* (J_serca + ks * RyRo * Caj * Vmyo / Vsr
//							+ (1.0 + 0.25 * AF) * 5.348e-6 * Caj
//							+ koff_csqn * Csqnb))
//			/ (1
//					+ dt
//							* (ks * RyRo * Vmyo / Vsr
//									+ (1.0 + 0.25 * AF) * 5.348e-6
//									+ kon_csqn * (Bmax_Csqn - Csqnb)));
	// ydot(31)=0;

	// Sodium Concentrations
	I_Na_tot_junc = I_Na_junc + I_nabk_junc + 3 * I_ncx_junc + 3 * I_nak_junc
			+ I_CaNa_junc + I_NaL_junc;   // [uA/uF]
	I_Na_tot_sl = I_Na_sl + I_nabk_sl + 3 * I_ncx_sl + 3 * I_nak_sl + I_CaNa_sl
			+ I_NaL_sl;   // [uA/uF]
	I_Na_tot_sl2 = 3 * I_ncx_sl + 3 * I_nak_sl + I_CaNa_sl;   // [uA/uF]
	I_Na_tot_junc2 = 3 * I_ncx_junc + 3 * I_nak_junc + I_CaNa_junc;   // [uA/uF]

	Naj += dt
			* (-I_Na_tot_junc * Cmem / (Vjunc * Frdy)
					+ J_na_juncsl / Vjunc * (Nasl - Naj) - dNaBj);
	Nasl += dt
			* (-I_Na_tot_sl * Cmem / (Vsl * Frdy)
					+ J_na_juncsl / Vsl * (Naj - Nasl)
					+ J_na_slmyo / Vsl * (Nai - Nasl) - dNaBsl);
	//FluxNaSL=ydot(33);
	// ydot(32) = 0;
	// ydot(33) = 0;
	Nai += dt * (J_na_slmyo / Vmyo * (Nasl - Nai));             // [mM/msec]
	// ydot(34)=0;

	// Potassium Concentration
	I_K_tot = I_to + I_kr + I_ks + I_ki - 2 * I_nak + I_CaK + I_kp + I_kur; // [uA/uF] //SVP: added IKur
	// ydot(35) = 0; //-I_K_tot*Cmem/(Vmyo*Frdy);           // [mM/msec]
	Ki += dt * 0.0;    //-I_K_tot*Cmem/(Vmyo*Frdy);

	// Calcium Concentrations
	I_Ca_tot_junc = I_Ca_junc + I_cabk_junc + I_pca_junc - 2 * I_ncx_junc; // [uA/uF]
	I_Ca_tot_sl = I_Ca_sl + I_cabk_sl + I_pca_sl - 2 * I_ncx_sl;      // [uA/uF]
	Caj += dt
			* (-I_Ca_tot_junc * Cmem / (Vjunc * 2 * Frdy)
					+ J_ca_juncsl / Vjunc * (Casl - Caj) - J_CaB_junction
					+ (J_SRCarel) * Vsr / Vjunc + J_SRleak * Vmyo / Vjunc); // Ca_j
	Casl += dt
			* (-I_Ca_tot_sl * Cmem / (Vsl * 2 * Frdy)
					+ J_ca_juncsl / Vsl * (Caj - Casl)
					+ J_ca_slmyo / Vsl * (Cai - Casl) - J_CaB_sl);   // Ca_sl
	// ydot(36)=0;
	// ydot(37)=0;
//    std::cout << "* Kharche11: Cai 1: " << Cai << std::flush;
	// ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol;//+J_ca_slmyo/Vmyo*(Casl-Cai);    // [mM/msec]
	Cai += dt
			* (-J_serca * Vsr / Vmyo - J_CaB_cytosol
					+ J_ca_slmyo / Vmyo * (Casl - Cai));
//    std::cout << ", Cai 2: " << Cai << std::endl;
	// ydot(38)=0;

	I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          // [uA/uF]
	I_Cl_tot = I_ClCa + I_Clbk + I_ClCFTR;                        // [uA/uF]
	I_Ca_tot = I_Ca_tot_junc + I_Ca_tot_sl;
	I_tot = I_Na_tot + I_Cl_tot + I_Ca_tot + I_K_tot;
}

} /* namespace BeatIt */
