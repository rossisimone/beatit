/*
 * NegroniLascano.cpp
 *
 *  Created on: Jul 13, 2017
 *      Author: srossi
 */
#include "Electromechanics/Activation/NegroniLascano.hpp"
#include <cmath>
namespace BeatIt
{

NegroniLascano::NegroniLascano() : super(5, "NegroniLascano")
{
    // Without potential
    M_variablesNames[0] = "TCa";
    M_variablesNames[1] = "TCas";
    M_variablesNames[2] = "Ts";
    M_variablesNames[3] = "X";
    M_variablesNames[4] = "F";

}

NegroniLascano::~NegroniLascano()
{
	// TODO Auto-generated destructor stub
}


void
NegroniLascano::solve(std::vector<double>& variables, double Cai, double dt, double I4f, double Cai_diast)
{

//	 Cai =  (*ionic_model_system.solution)(var_index);
//	 Cai /= Cai_rescale;
//	 if(Cai <= Cai_diastolic) Cai = 0.0;
//	 else Cai -= Cai_diastolic;
	 TCa  = variables[0];
	 TCas = variables[1];
	 Ts   = variables[2];
	 X    = variables[3];
	 F    = variables[4];

	 TCa_eff  = TCa * std::exp(-R*(I4f - La)*(I4f - La)) ;
	 T = Tt - TCa - TCas - Ts;
	 dX = B * (I4f - X - hc); // dX/dt
	 Qb = Y1 * Cai * T - Z1 * TCa;
	 Qa = Y2 * TCa_eff - Z2 * TCas;
	 Qr = Y3 * TCas - Z3 * Ts * Cai ;
	 Qd = Y4 * Ts ;
	 Qd1 = Yd * dX * dX * Ts ;
	 Qd2 = Yd * dX * dX * TCas ;
	 dTCa = Qb - Qa ; // dTCa/dt
	 dTCas = Qa - Qr - Qd2 ; // dTCa*/dt
	 dTs = Qr - Qd - Qd1 ; // dT*/dt
	 TCa +=  dt * dTCa;
	 TCas+= dt * dTCas;
	 Ts += dt * dTs;
	 X += dt * dX;

     F = A * ( TCas + Ts) * (I4f-X);

     variables[0] = TCa;
	 variables[1] = TCas;
	 variables[2] = Ts;
	 variables[3] = X;
	 variables[4] = F;
}
//! Initialize the values of the variables
/*!
 *  \param [in] variables Vector containing the local value of all variables
 */
void
NegroniLascano::initialize(std::vector<double>& variables)
{

}

//! Initialize the output files with the names
/*!
 *  \param [in] output output file where we will save the values
 */
void
NegroniLascano::initializeSaveData(std::ostream& output)
{
    // time -  0
    output << "time TCa TCas Ts X F\n";
}




} /* namespace BeatIt */
