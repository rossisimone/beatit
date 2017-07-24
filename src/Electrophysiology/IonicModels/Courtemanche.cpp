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
 * \file Courtemanche.cpp
 *
 * \class Courtemanche
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

#include "Electrophysiology/IonicModels/Courtemanche.hpp"
#include <cmath>

namespace BeatIt
{

IonicModel* createCourtemanche()
{
	return new Courtemanche;
}

Courtemanche::Courtemanche() :
		super(27, 0, "Courtemanche", CellType::MCell)
{
	// The name refer to the comments below the lines
	// Nai = 11.2
	M_variablesNames[0] = "Nai";
	//Ki = 139
	M_variablesNames[1] = "Ki";
	//Cai = 0.000102
	M_variablesNames[2] = "Cai";
	// m0 = 0.00291
	M_variablesNames[3] = "m";
	//ho=0.965;
	M_variablesNames[4] = "h";
	// j = 0.978;
	M_variablesNames[5] = "j";
	//d = 0.000137;
	M_variablesNames[6] = "d";
	// f = 0.999837;
	M_variablesNames[7] = "f";
	//xs = 0.0187;
	M_variablesNames[8] = "xs";
	//xr = 0.0000329;
	M_variablesNames[9] = "xr";
	//ato = 0.0304;
	M_variablesNames[10] = "ato";
	//iito = 0.999;
	M_variablesNames[11] = "iito";
	//uakur = 0.00496;
	M_variablesNames[12] = "uakur";
	//uikur = 0.999;
	M_variablesNames[13] = "uikur";
	//fca = 0.775;
	M_variablesNames[14] = "fca";
	//ireljsrol=0;
	M_variablesNames[15] = "ireljsrol";
	//jsr = 1.49;
	M_variablesNames[16] = "jsr";
	//nsr = 1.49;
	M_variablesNames[17] = "nsr";
	//trpn = 0.0118;
	M_variablesNames[18] = "trpn";
	//cmdn = 0.00205;
	M_variablesNames[19] = "cmdn";
	//csqn = 6.51;
	M_variablesNames[20] = "csqn";
	//urel = 0.00;
	M_variablesNames[21] = "urel";
	//vrel = 1.00;
	M_variablesNames[22] = "vrel";
	//wrel = 0.999;
	M_variablesNames[23] = "wrel";
	//yach = 2.54e-2;
	M_variablesNames[24] = "yach";
	//iky = 0.6;
	M_variablesNames[25] = "iky";

}

void Courtemanche::initializeSaveData(std::ostream& output)
{
	// time -  0
	output << "time v ";
	for (auto && var_name : M_variablesNames)
	{
		output << var_name << " ";
	}
	output << "\n";
}

void Courtemanche::initialize(std::vector<double>& variables)
{
	// v = -81.2;       /* Initial Voltage (mv) */
	variables[0] = -81.2;
	// Nai = 11.2(mM)
	variables[1] = 11.2;
	//Ki = 139(mM)
	variables[2] = 139;
	//Cai = 0.000102 (mM)
	variables[3] = 0.000102;
	// m0 = 0.00291
	variables[4] = 0.00291;
	//ho=0.965;
	variables[5] = 0.965;
	// j = 0.978;
	variables[6] = 0.978;
	//d = 0.000137;
	variables[7] = 0.000137;
	// f = 0.999837;
	variables[8] = 0.999837;
	//xs = 0.0187;
	variables[9] = 0.0187;
	//xr = 0.0000329;
	variables[10] = 0.0000329;
	//ato = 0.0304;
	variables[11] = 0.0304;
	//iito = 0.999;
	variables[12] = 0.999;
	//uakur = 0.00496;
	variables[13] = 0.00496;
	//uikur = 0.999;
	variables[14] = 0.999;
	//fca = 0.775;
	variables[15] = 0.775;
	//ireljsrol=0;
	variables[16] = 0;
	//jsr = 1.49;
	variables[17] = 1.49;
	//nsr = 1.49;
	variables[18] = 1.49;
	//trpn = 0.0118;
	variables[19] = 0.0118;
	//cmdn = 0.00205;
	variables[20] = 0.00205;
	//csqn = 6.51;
	variables[21] = 6.51;
	//urel = 0.00;
	variables[22] = 0.0;
	//vrel = 1.00;
	variables[23] = 1.0;
	//wrel = 0.999;
	variables[24] = 0.999;
	//yach = 2.54e-2;
	variables[25] = 2.54e-2;
	//iky = 0.6;
	variables[26] = 0.6;
}

//! Solve method
/*!
 *  \param [in] variables Vector containing the local value of all variables
 *  \param [in] appliedCurrent value of the applied current
 *  \param [in] dt        Timestep
 */
void Courtemanche::solve(std::vector<double>& variables, double appliedCurrent,
		double dt)
{
	updateVariables(variables, appliedCurrent, dt);
	variables[0] += dt * evaluateIonicCurrent(variables, appliedCurrent, dt);
}

//! Evaluate total ionic current for the computation of the potential
/*!
 *  \param [in] variables Vector containing the local value of all variables
 *  \param [in] appliedCurrent value of the applied current
 *  \param [in] dt        Timestep
 */
double Courtemanche::evaluateIonicCurrent(std::vector<double>& variables,
		double appliedCurrent, double dt)
{
	return -it;
}

void Courtemanche::updateVariables(std::vector<double>& variables,
		double appliedCurrent, double dt)
{
	this->dt = dt;
	st = -appliedCurrent;
	v = variables[0];
	/*  Ion Concentrations */
	nai = variables[1];
	ki = variables[2];
	cai = variables[3];

	/*  Gate Conditions */
	m = variables[4];
	h = variables[5];
	j = variables[6];
	d = variables[7];
	f = variables[8];
	xs = variables[9];
	xr = variables[10];
	ato = variables[11];
	iito = variables[12];
	uakur = variables[13];
	uikur = variables[14];
	fca = variables[15];
	ireljsrol = variables[16];

	jsr = variables[17];
	nsr = variables[18];
	trpn = variables[19];
	cmdn = variables[20];
	csqn = variables[21];
	urel = variables[22];
	vrel = variables[23];
	wrel = variables[24];
	yach = variables[25];
	iky = variables[26];

    update();
	/*  Ion Concentrations */
	variables[1] = nai;
	variables[2] = ki;
	variables[3] = cai;

	/*  Gate Conditions */
	variables[4] = m;
	variables[5] = h;
	variables[6] = j;
	variables[7] = d;
	variables[8] = f;
	variables[9] = xs;
	variables[10] = xr;
	variables[11] = ato;
	variables[12] = iito;
	variables[13] = uakur;
	variables[14] = uikur;
	variables[15] = fca;
	variables[16] = ireljsrol;

	variables[17] = jsr;
	variables[18] = nsr;
	variables[19] = trpn;
	variables[20] = cmdn;
	variables[21] = csqn;
	variables[22] = urel;
	variables[23] = vrel;
	variables[24] = wrel;
	variables[25] = yach;
	variables[26] = iky;
}


/********************************************************/
void
Courtemanche::update()
{
	            comp_ina ();
                comp_ical ();
                comp_ikr ();
                comp_iks ();
                comp_iki ();
				comp_ikach ();
                comp_ikur ();
				comp_ito ();
                comp_inaca ();
                comp_inak ();
                comp_ipca ();
                comp_icab ();
                comp_inab ();
                comp_it ();

                conc_nai ();
                conc_ki ();
                conc_nsr ();
                conc_jsr ();
                calc_itr ();
                conc_cai ();
}

/* Functions that describe the currents begin here */

void Courtemanche::comp_ina()
{
	gna = 7.8;
	ena = ((R * temp) / frdy) * log(nao / nai);

	am = 0.32 * (v + 47.13) / (1 - exp(-0.1 * (v + 47.13)));
	bm = 0.08 * exp(-v / 11);

	if (v < -40)
	{
		ah = 0.135 * exp((80 + v) / -6.8);
		bh = 3.56 * exp(0.079 * v) + 310000 * exp(0.35 * v);
		aj = (-127140 * exp(0.2444 * v) - 0.00003474 * exp(-0.04391 * v))
				* ((v + 37.78) / (1 + exp(0.311 * (v + 79.23))));
		bj = (0.1212 * exp(-0.01052 * v)) / (1 + exp(-0.1378 * (v + 40.14)));
	}

	else
	{
		ah = 0.0;
		bh = 1 / (0.13 * (1 + exp((v + 10.66) / -11.1)));
		aj = 0.0;
		bj = (0.3 * exp(-0.0000002535 * v)) / (1 + exp(-0.1 * (v + 32)));
	}

	h = ah / (ah + bh) - ((ah / (ah + bh)) - h) * exp(-dt / (1 / (ah + bh)));
	j = aj / (aj + bj) - ((aj / (aj + bj)) - j) * exp(-dt / (1 / (aj + bj)));
	m = am / (am + bm) - ((am / (am + bm)) - m) * exp(-dt / (1 / (am + bm)));

	ina = gna * m * m * m * h * j * (v - ena);
}

void Courtemanche::comp_ical()
{
	dss = 1 / (1 + exp(-(v + 10) / 8));
	taud = (1 - exp((v + 10) / -6.24))
			/ (0.035 * (v + 10) * (1 + exp((v + 10) / -6.24)));

	fss = 1 / (1 + exp((v + 28) / 6.9));
	tauf = 9 / (0.0197 * exp(-pow((0.0337 * (v + 10)), 2)) + 0.02);

	fcass = 1 / (1 + cai / 0.00035);
	taufca = 2;

	d = dss - (dss - d) * exp(-dt / taud);
	f = fss - (fss - f) * exp(-dt / tauf);
	fca = fcass - (fcass - fca) * exp(-dt / tauf);

	ibarca = gcalbar * (v - 65);

	ilca = d * f * fca * ibarca;

	ilcatot = ilca;
}

void Courtemanche::comp_ikr()
{
	gkr = 0.0294 * sqrt(ko / 5.4);
	ekr = ((R * temp) / frdy) * log(ko / ki);

	xrss = 1 / (1 + exp(-(v + 14.1) / 6.5));
	tauxr = 1
			/ (0.0003 * (v + 14.1) / (1 - exp(-(v + 14.1) / 5))
					+ 0.000073898 * (v - 3.3328)
							/ (exp((v - 3.3328) / 5.1237) - 1));

	xr = xrss - (xrss - xr) * exp(-dt / tauxr);

	r = 1 / (1 + exp((v + 15) / 22.4));

	ikr = gkr * xr * r * (v - ekr);
}

void Courtemanche::comp_iks()
{
	gks = 0.129;
	eks = ((R * temp) / frdy) * log(ko / ki);
	tauxs = 0.5
			/ (0.00004 * (v - 19.9) / (1 - exp(-(v - 19.9) / 17))
					+ 0.000035 * (v - 19.9) / (exp((v - 19.9) / 9) - 1));
	xsss = 1 / pow((1 + exp(-(v - 19.9) / 12.7)), 0.5);
	xs = xsss - (xsss - xs) * exp(-dt / tauxs);

	iks = gks * xs * xs * (v - eks);
}

void Courtemanche::comp_iki()
{
	gki = 0.09 * pow(ko / 5.4, 0.4);
	eki = ((R * temp) / frdy) * log(ko / ki);

	kin = 1 / (1 + exp(0.07 * (v + 80)));

	iki = gki * kin * (v - eki);
	/*

	 // modified from Matsuoka, et al Jap J Physiol 2003:53:105-123
	 iku = 0.75*exp(0.035*(v-eki-10))/(1+exp(0.015*(v-eki-140)));
	 ikl = 3*exp(-0.048*(v-eki-10))*(1+exp(0.064*(v-eki-38)))/(1+exp(0.03*(v-eki-70)));
	 ikay =1/(8000*exp((v-eki-97)/8.5)+7*exp((v-eki-97)/300));
	 ikby =1/(0.00014*exp(-(v-eki-97)/9.1)+0.2*exp(-(v-eki-97)/500));
	 tauiky = 1/(ikay+ikby);
	 ikyss = ikay/(ikay+ikby);
	 iky = ikyss - (ikyss-iky)*exp(-dt/tauiky);
	 foiki = ikl/(iku+ikl);
	 fbiki = iku/(iku+ikl);


	 iki = gki*(pow(foiki,4)+8*pow(foiki,3)*fbiki/3+2*foiki*foiki*fbiki*fbiki)*iky*(v-eki);
	 */
}

void Courtemanche::comp_ikach()
{
	gkach = 0.135;
	ekach = ((R * temp) / frdy) * log(ko / ki);
	alphayach = 1.232e-2 / (1 + 0.0042 / ach) + 0.0002475;
	betayach = 0.01 * exp(0.0133 * (v + 40));
	tauyach = 1 / (alphayach + betayach);
	yachss = alphayach / (alphayach + betayach);

	yach = yachss - (yachss - yach) * exp(-dt / tauyach);
	ikach = gkach * yach * (v - ekach) / (1 + exp((v + 20) / 20));
}

void Courtemanche::comp_ikur()
{
	gkur = 0.005 + 0.05 / (1 + exp(-(v - 15) / 13));
	ekur = ((R * temp) / frdy) * log(ko / ki);
	alphauakur = 0.65 / (exp(-(v + 10) / 8.5) + exp(-(v - 30) / 59.0));
	betauakur = 0.65 / (2.5 + exp((v + 82) / 17.0));
	tauuakur = 1 / (3 * (alphauakur + betauakur));
	uakurss = 1 / (1 + exp(-(v + 30.3) / 9.6));
	alphauikur = 1 / (21 + exp(-(v - 185) / 28));
	betauikur = exp((v - 158) / 16);
	tauuikur = 1 / (3 * (alphauikur + betauikur));
	uikurss = 1 / (1 + exp((v - 99.45) / 27.48));

	uakur = uakurss - (uakurss - uakur) * exp(-dt / tauuakur);
	uikur = uikurss - (uikurss - uikur) * exp(-dt / tauuikur);

	ikur = gkur * uakur * uakur * uakur * uikur * (v - ekur);
}

void Courtemanche::comp_ito()
{
	gito = 0.1652;
	erevto = ((R * temp) / frdy) * log(ko / ki);

	alphaato = 0.65 / (exp(-(v + 10) / 8.5) + exp(-(v - 30) / 59));
	betaato = 0.65 / (2.5 + exp((v + 82) / 17));
	tauato = 1 / (3 * (alphaato + betaato));
	atoss = 1 / (1 + exp(-(v + 20.47) / 17.54));
	ato = atoss - (atoss - ato) * exp(-dt / tauato);

	alphaiito = 1 / (18.53 + exp((v + 113.7) / 10.95));
	betaiito = 1 / (35.56 + exp(-(v + 1.26) / 7.44));
	tauiito = 1 / (3 * (alphaiito + betaiito));
	iitoss = 1 / (1 + exp((v + 43.1) / 5.3));
	iito = iitoss - (iitoss - iito) * exp(-dt / tauiito);

	ito = gito * ato * ato * ato * iito * (v - erevto);
}

void Courtemanche::comp_inaca()
{
	inaca =
			1750
					* (exp(gammas * frdy * v / (R * temp)) * nai * nai * nai
							* cao
							- exp((gammas - 1) * frdy * v / (R * temp)) * nao
									* nao * nao * cai)
					/ ((pow(kmnancx, 3) + pow(nao, 3)) * (kmcancx + cao)
							* (1
									+ ksatncx
											* exp(
													(gammas - 1) * frdy * v
															/ (R * temp))));
}

void Courtemanche::comp_inak()
{
	sigma = (exp(nao / 67.3) - 1) / 7;

	//fnak = 1/(1+0.1245*exp((-0.1*v*frdy)/(R*temp))+0.0365*sigma*exp((-v*frdy)/(R*temp)));
	fnak = (v + 150) / (v + 200);
	inak = ibarnak * fnak * (1 / (1 + pow((kmnai / nai), 1.5)))
			* (ko / (ko + kmko));
}

void Courtemanche::comp_ipca()
{
	ipca = (ibarpca * cai) / (kmpca + cai);
}

void Courtemanche::comp_icab()
{
	gcab = 0.00113;
	ecan = ((R * temp) / frdy) * log(cao / cai);

	icab = gcab * (v - ecan);
}

void Courtemanche::comp_inab()
{
	gnab = 0.000674;
	enan = ((R * temp) / frdy) * log(nao / nai);

	inab = gnab * (v - enan);
}

/* Total sum of currents is calculated here, if the time is between stimtime = 0 and stimtime = 0.5, a stimulus is applied */
void Courtemanche::comp_it()
{
	naiont = ina + inab + 3 * inak + 3 * inaca + 1.5e-2;
	kiont = ikr + iks + iki - 2 * inak + ito + ikur + ikach + 1.5e-2;
	caiont = ilca + icab + ipca - 2 * inaca;
	it = st + naiont + kiont + caiont;
}

/* Functions that calculate intracellular ion concentrations begins here */

void Courtemanche::conc_nai()
{
	dnai = -dt * naiont * acap / (vmyo * zna * frdy);
	nai = dnai + nai;
}

void Courtemanche::conc_ki()
{
	dki = -dt * kiont * acap / (vmyo * zk * frdy);
	ki = dki + ki;
}

void Courtemanche::conc_nsr()
{
	kleak = iupbar / nsrbar;
	ileak = kleak * nsr;

	iup = iupbar * cai / (cai + kmup);
	csqn = csqnbar * (jsr / (jsr + kmcsqn));

	dnsr = dt * (iup - ileak - itr * vjsr / vnsr);
	nsr = dnsr + nsr;
}

void Courtemanche::conc_jsr()
{

//		fn = vjsr*(1e-12)*ireljsrol-(5e-13)*(ilca/2+inaca/5)*acap/frdy;
	fn = vjsr * (1e-12) * ireljsrol - (1e-12) * caiont * acap / (2 * frdy);

	tauurel = 8.0;
	urelss = 1 / (1 + exp(-(fn - 3.4175e-13) / 13.67e-16));
	tauvrel = 1.91 + 2.09 / (1 + exp(-(fn - 3.4175e-13) / 13.67e-16));
	vrelss = 1 - 1 / (1 + exp(-(fn - 6.835e-14) / 13.67e-16));
	tauwrel = 6.0 * (1 - exp(-(v - 7.9) / 5))
			/ ((1 + 0.3 * exp(-(v - 7.9) / 5)) * (v - 7.9));
	wrelss = 1 - 1 / (1 + exp(-(v - 40) / 17));

	urel = urelss - (urelss - urel) * exp(-dt / tauurel);
	vrel = vrelss - (vrelss - vrel) * exp(-dt / tauvrel);
	wrel = wrelss - (wrelss - wrel) * exp(-dt / tauwrel);

	greljsrol = grelbarjsrol * urel * urel * vrel * wrel;
	ireljsrol = greljsrol * (jsr - cai);

	djsr = dt * (itr - 0.5 * ireljsrol)
			/ (1 + csqnbar * kmcsqn / pow((jsr + kmcsqn), 2)); //LAI

	jsr = djsr + jsr;
}

void Courtemanche::calc_itr()
{
	itr = (nsr - jsr) / tautr;
}

void Courtemanche::conc_cai()
{
	trpn = trpnbar * (cai / (cai + kmtrpn));
	cmdn = cmdnbar * (cai / (cai + kmcmdn));

	b1cai = -caiont * acap / (2 * frdy * vmyo)
			+ (vnsr * (ileak - iup) + 0.5 * ireljsrol * vjsr) / vmyo; //LAI
	b2cai = 1 + trpnbar * kmtrpn / pow((cai + kmtrpn), 2)
			+ cmdn * kmcmdn / pow((cai + kmcmdn), 2);
	dcai = dt * b1cai / b2cai;

	cai = dcai + cai;
}

} /* namespace BeatIt */
