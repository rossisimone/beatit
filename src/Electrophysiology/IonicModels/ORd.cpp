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
 * \file ORd.cpp
 *
 * \class ORd
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
 * Created on: Aug 7, 2016
 *
 */

#include <cmath>
#include <fstream>
#include "Electrophysiology/IonicModels/ORd.hpp"

namespace BeatIt
{


IonicModel* createORd()
{
	return new ORd();
}

ORd::ORd()
  : super(41,0, "ORd", CellType::MCell)
{
    // Without potential

    //nai=7;
    M_variablesNames[0] = "nai";
    //nass=nai;
    M_variablesNames[1] = "nass";
    //ki=145;
    M_variablesNames[2] = "ki";
    //kss=ki;
    M_variablesNames[3] = "kss";
    //cai=1.0e-4;
    M_variablesNames[4] = "cai";
    //cass=cai;
    M_variablesNames[5] = "cass";
    //cansr=1.2;
    M_variablesNames[6] = "cansr";
    //cajsr=cansr;
    M_variablesNames[7] = "caisr";
    //m=0;
    M_variablesNames[8] = "m";
    //hf=1;
    M_variablesNames[9] = "hf";
    //hs=1;
    M_variablesNames[10] = "hs";
    //j=1;
    M_variablesNames[11] = "j";
    //hsp=1;
    M_variablesNames[12] = "hsp";
    //jp=1;
    M_variablesNames[13] = "jp";
    //mL=0;
    M_variablesNames[14] = "mL";
    //hL=1;
    M_variablesNames[15] = "hL";
    //hLp=1;
    M_variablesNames[16] = "hLp";
    //a=0;
    M_variablesNames[17] = "a";
    //iF=1;
    M_variablesNames[18] = "iF";
    //iS=1;
    M_variablesNames[19] = "iS";
    //ap=0;
    M_variablesNames[20] = "ap";
    //iFp=1;
    M_variablesNames[21] = "iFp";
    //iSp=1;
    M_variablesNames[22] = "iSp";
    //d=0;
    M_variablesNames[23] = "d";
    //ff=1;
    M_variablesNames[24] = "ff";
    //fs=1;
    M_variablesNames[25] = "fs";
    //fcaf=1;
    M_variablesNames[26] = "fcaf";
    //fcas=1;
    M_variablesNames[27] = "fcas";
    //jca=1;
    M_variablesNames[28] = "jca";
    //nca=0;
    M_variablesNames[29] = "nca";
    //ffp=1;
    M_variablesNames[30] = "ffp";
    //fcafp=1;
    M_variablesNames[31] = "fcafp";
    //xrf=0;
    M_variablesNames[32] = "xrf";
    //xrs=0;
    M_variablesNames[33] = "xrs";
    //xs1=0;
    M_variablesNames[34] = "xs1";
    //xs2=0;
    M_variablesNames[35] = "xs2";
    //xk1=1;
    M_variablesNames[36] = "xk1";
    //Jrelnp=0;
    M_variablesNames[37] = "Jrelnp";
    //Jrelp=0;
    M_variablesNames[38] = "Jrelp";
    //CaMKt=0;
    M_variablesNames[39] = "CaMKt";
}

void
ORd::initialize(std::vector<double>& variables)
{
    // V
	variables[0] = -87.5;
    //nai=7;
    variables[1] = 7.;
    //nass=nai;
    variables[2] = 7.;
    //ki=145;
    variables[3] = 145.;
    //kss=ki;
    variables[4] = 145.;
    //cai=1.0e-4;
    variables[5] = 1e-4;
    //cass=cai;
    variables[6] = 1e-4;
    //cansr=1.2;
    variables[7] = 1.2;
    //cajsr=cansr;
    variables[8] = 1.2;
    //m=0;
    variables[9] = 0;
    //hf=1;
    variables[10] = 1;
    //hs=1;
    variables[11] = 1;
    //j=1;
    variables[12] = 1;
    //hsp=1;
    variables[13] = 1;
    //jp=1;
    variables[14] = 1;
    //mL=0;
    variables[15] = 0;
    //hL=1;
    variables[16] = 1;
    //hLp=1;
    variables[17] = 1;
    //a=0;
    variables[18] = 0;
    //iF=1;
    variables[19] = 1;
    //iS=1;
    variables[20] = 1;
    //ap=0;
    variables[21] = 0;
    //iFp=1;
    variables[22] = 1;
    //iSp=1;
    variables[23] = 1;
    //d=0;
    variables[24] = 0;
    //ff=1;
    variables[25] = 1;
    //fs=1;
    variables[26] = 1;
    //fcaf=1;
    variables[27] = 1;
    //fcas=1;
    variables[28] = 1;
    //jca=1;
    variables[29] = 1;
    //nca=0;
    variables[30] = 0;
    //ffp=1;
    variables[31] = 1;
    //fcafp=1;
    variables[32] = 1;
    //xrf=0;
    variables[33] = 0;
    //xrs=0;
    variables[34] = 0;
    //xs1=0;
    variables[35] = 0;
    //xs2=0;
    variables[36] = 0;
    //xk1=1;
    variables[37] = 1;
    //Jrelnp=0;
    variables[38] = 0;
    //Jrelp=0;
    variables[39] = 0;
    //CaMKt=0;
    variables[40] = 0;

}


//! Solve method
/*!
 *  \param [in] variables Vector containing the local value of all variables
 *  \param [in] appliedCurrent value of the applied current
 *  \param [in] dt        Timestep
 */
void
ORd::solve(std::vector<double>& variables, double appliedCurrent, double dt)
{
	updateVariables(variables, appliedCurrent, dt);
    variables[0] += dt * evaluateIonicCurrent(variables, appliedCurrent, dt);
}

//! Update all the variables in the ionic model
/*!
 *  \param [in] variables Vector containing the local value of all variables
 *  \param [in] dt        Timestep
 */
void
ORd::updateVariables(std::vector<double>& variables, double appliedCurrent, double dt)
{
	// For compatibility  with the original code where the applied stimulus in opposite
	Ist = -appliedCurrent;
    revpots(variables);
    RGC(variables, dt);
    FBC(variables, dt);
}

//! Evaluate total ionic current for the computation of the potential
/*!
 *  \param [in] variables Vector containing the local value of all variables
 *  \param [in] appliedCurrent value of the applied current
 *  \param [in] dt        Timestep
 */
double
ORd::evaluateIonicCurrent(std::vector<double>& variables, double appliedCurrent, double dt)
{
	// For compatibility  with the original code where the applied stimulus in opposite
	Ist =-appliedCurrent;
	return -(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa+INaK+INab+IKb+IpCa+ICab+Ist);

}


void
ORd::revpots(std::vector<double>& variables)
{
    double nai = variables[1];
    double ki  = variables[3];
    ENa=(R*T/F)*std::log(nao/nai);
    EK=(R*T/F)*std::log( ko/ki);
    EKs=(R*T/F)*std::log(( ko+0.01833* nao)/(ki+0.01833*nai));
}

void
ORd::initializeSaveData(std::ostream& output)
{
	// time -  0
	output << "time v ";
	//  1 - 10
	output << "nai nass ki kss cai cass cansr cajsr m hf ";
	// 11 - 20
	output << "hs j hsp jp mL hL hLp a iF iS ";
	// 21 - 30
	output << "ap iFp iSp d ff fs fcaf fcas jca nca";
	// 31 - 40
	output << "ffp fcafp xrf xrs xs1 xs2 xk1 Jrelnp Jrelp CaMKt\n";
}

void
ORd::RGC(std::vector<double>& variables, double dt)
{
	double v     = variables[0];
	double nai   = variables[1];
	double nass  = variables[2];
	double ki    = variables[3];
	double kss   = variables[4];
	double cai   = variables[5];
	double cass  = variables[6];
	double cansr = variables[7];
	double cajsr = variables[8];
	double m     = variables[9];
	double hf    = variables[10];
	double hs    = variables[11];
	double j     = variables[12];
	double hsp   = variables[13];
	double jp    = variables[14];
	double mL    = variables[15];
	double hL    = variables[16];
	double hLp   = variables[17];
	double a     = variables[18];
	double iF    = variables[19];
	double iS    = variables[20];
	double ap    = variables[21];
	double iFp   = variables[22];
	double iSp   = variables[23];
	double d     = variables[24];
	double ff    = variables[25];
	double fs    = variables[26];
	double fcaf  = variables[27];
	double fcas  = variables[28];
	double jca   = variables[29];
	double nca   = variables[30];
	double ffp   = variables[31];
	double fcafp = variables[32];
	double xrf   = variables[33];
	double xrs   = variables[34];
	double xs1   = variables[35];
	double xs2   = variables[36];
	double xk1   = variables[37];
	//double Jrelnp= variables[38];
	//double Jrelp = variables[39];
	double CaMKt = variables[40];


	CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
	CaMKa=CaMKb+CaMKt;
	double vffrt=v*F*F/(R*T);
	double vfrt=v*F/(R*T);

	// m variable
	double mss=1.0/(1.0+std::exp((-(v+39.57))/9.871));
	double tm=1.0/(6.765*std::exp((v+11.64)/34.77)+8.552*std::exp(-(v+77.42)/5.955));
	m=mss-(mss-m)*std::exp(-dt/tm);

	// h variable
	double hss=1.0/(1+std::exp((v+82.90)/6.086));
	double thf=1.0/(1.432e-5*std::exp(-(v+1.196)/6.285)+6.149*std::exp((v+0.5096)/20.27));
	double ths=1.0/(0.009794*std::exp(-(v+17.95)/28.05)+0.3343*std::exp((v+5.730)/56.66));
	double Ahf=0.99;
	double Ahs=1.0-Ahf;
	hf=hss-(hss-hf)*std::exp(-dt/thf);
	hs=hss-(hss-hs)*std::exp(-dt/ths);

	double h=Ahf*hf+Ahs*hs;
	double jss=hss;
	double tj=2.038+1.0/(0.02136*std::exp(-(v+100.6)/8.281)+0.3052*std::exp((v+0.9941)/38.45));
	j=jss-(jss-j)*std::exp(-dt/tj);

	double hssp=1.0/(1+std::exp((v+89.1)/6.086));
	double thsp=3.0*ths;
	hsp=hssp-(hssp-hsp)*std::exp(-dt/thsp);

	double hp=Ahf*hf+Ahs*hsp;
	double tjp=1.46*tj;
	jp=jss-(jss-jp)*std::exp(-dt/tjp);

	double GNa=75;
	double fINap=(1.0/(1.0+KmCaMK/CaMKa));
	INa=GNa*(v-ENa)*m*m*m*((1.0-fINap)*h*j+fINap*hp*jp);

	double mLss=1.0/(1.0+std::exp((-(v+42.85))/5.264));
	double tmL=tm;
	mL=mLss-(mLss-mL)*std::exp(-dt/tmL);

	double hLss=1.0/(1.0+std::exp((v+87.61)/7.488));
	double thL=200.0;
	hL=hLss-(hLss-hL)*std::exp(-dt/thL);

	double hLssp=1.0/(1.0+std::exp((v+93.81)/7.488));
	double thLp=3.0*thL;
	hLp=hLssp-(hLssp-hLp)*std::exp(-dt/thLp);

	double GNaL=0.0075;
	if (M_cellType == CellType::Epicardial )
	{
		GNaL*=0.6;
	}
	double fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
	INaL=GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);

	double ass=1.0/(1.0+std::exp((-(v-14.34))/14.82));
	double ta=1.0515/(1.0/(1.2089*(1.0+std::exp(-(v-18.4099)/29.3814)))+3.5/(1.0+std::exp((v+100.0)/29.3814)));
	a=ass-(ass-a)*std::exp(-dt/ta);

	double iss=1.0/(1.0+std::exp((v+43.94)/5.711));
	double delta_epi;
	if (M_cellType == CellType::Epicardial)
	{
		delta_epi=1.0-(0.95/(1.0+std::exp((v+70.0)/5.0)));
	}
	else
	{
		delta_epi=1.0;
	}

	double tiF=4.562+1/(0.3933*std::exp((-(v+100.0))/100.0)+0.08004*std::exp((v+50.0)/16.59));
	double tiS=23.62+1/(0.001416*std::exp((-(v+96.52))/59.05)+1.780e-8*std::exp((v+114.1)/8.079));
	tiF*=delta_epi;
	tiS*=delta_epi;

	double AiF=1.0/(1.0+std::exp((v-213.6)/151.2));
	double AiS=1.0-AiF;
	iF=iss-(iss-iF)*std::exp(-dt/tiF);
	iS=iss-(iss-iS)*std::exp(-dt/tiS);

	double i=AiF*iF+AiS*iS;
	double assp=1.0/(1.0+std::exp((-(v-24.34))/14.82));
	ap=assp-(assp-ap)*std::exp(-dt/ta);

	double dti_develop=1.354+1.0e-4/(std::exp((v-167.4)/15.89)+std::exp(-(v-12.23)/0.2154));
	double dti_recover=1.0-0.5/(1.0+std::exp((v+70.0)/20.0));
	double tiFp=dti_develop*dti_recover*tiF;
	double tiSp=dti_develop*dti_recover*tiS;
	iFp=iss-(iss-iFp)*std::exp(-dt/tiFp);
	iSp=iss-(iss-iSp)*std::exp(-dt/tiSp);

	double ip=AiF*iFp+AiS*iSp;
	double Gto=0.02;
	if (M_cellType == CellType::Epicardial)
	{
		Gto*=4.0;
	}
	else if (M_cellType == CellType::MCell )
	{
		Gto*=4.0;
	}
	double fItop=(1.0/(1.0+KmCaMK/CaMKa));
	Ito=Gto*(v-EK)*((1.0-fItop)*a*i+fItop*ap*ip);

	double dss=1.0/(1.0+std::exp((-(v+3.940))/4.230));
	double td=0.6+1.0/(std::exp(-0.05*(v+6.0))+std::exp(0.09*(v+14.0)));
	d=dss-(dss-d)*std::exp(-dt/td);

	double fss=1.0/(1.0+std::exp((v+19.58)/3.696));
	double tff=7.0+1.0/(0.0045*std::exp(-(v+20.0)/10.0)+0.0045*std::exp((v+20.0)/10.0));
	double tfs=1000.0+1.0/(0.000035*std::exp(-(v+5.0)/4.0)+0.000035*std::exp((v+5.0)/6.0));
	double Aff=0.6;
	double Afs=1.0-Aff;
	ff=fss-(fss-ff)*std::exp(-dt/tff);
	fs=fss-(fss-fs)*std::exp(-dt/tfs);

	double f=Aff*ff+Afs*fs;
	double fcass=fss;
	double tfcaf=7.0+1.0/(0.04*std::exp(-(v-4.0)/7.0)+0.04*std::exp((v-4.0)/7.0));
	double tfcas=100.0+1.0/(0.00012*std::exp(-v/3.0)+0.00012*std::exp(v/7.0));
	double Afcaf=0.3+0.6/(1.0+std::exp((v-10.0)/10.0));
	double Afcas=1.0-Afcaf;
	fcaf=fcass-(fcass-fcaf)*std::exp(-dt/tfcaf);
	fcas=fcass-(fcass-fcas)*std::exp(-dt/tfcas);

	double fca=Afcaf*fcaf+Afcas*fcas;
	double tjca=75.0;
	jca=fcass-(fcass-jca)*std::exp(-dt/tjca);

	double tffp=2.5*tff;
	ffp=fss-(fss-ffp)*std::exp(-dt/tffp);

	double fp=Aff*ffp+Afs*fs;
	double tfcafp=2.5*tfcaf;
	fcafp=fcass-(fcass-fcafp)*std::exp(-dt/tfcafp);

	double fcap=Afcaf*fcafp+Afcas*fcas;
	double Kmn=0.002;
	double k2n=1000.0;
	double km2n=jca*1.0;
	double anca=1.0/(k2n/km2n+std::pow(1.0+Kmn/cass,4.0));
	nca=anca*k2n/km2n-(anca*k2n/km2n-nca)*std::exp(-km2n*dt);

	double PhiCaL=4.0*vffrt*(cass*std::exp(2.0*vfrt)-0.341*cao)/(std::exp(2.0*vfrt)-1.0);
	double PhiCaNa=1.0*vffrt*(0.75*nass*std::exp(1.0*vfrt)-0.75*nao)/(std::exp(1.0*vfrt)-1.0);
	double PhiCaK=1.0*vffrt*(0.75*kss*std::exp(1.0*vfrt)-0.75*ko)/(std::exp(1.0*vfrt)-1.0);
	double zca=2.0;
	double PCa=0.0001;
	if (M_cellType == CellType::Epicardial)
	{
		PCa*=1.2;
	}
	else if (M_cellType == CellType::MCell)
	{
		PCa*=2.5;
	}
	double PCap=1.1*PCa;
	double PCaNa=0.00125*PCa;
	double PCaK=3.574e-4*PCa;
	double PCaNap=0.00125*PCap;
	double PCaKp=3.574e-4*PCap;
	double fICaLp=(1.0/(1.0+KmCaMK/CaMKa));
	ICaL=(1.0-fICaLp)*PCa*PhiCaL*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCap*PhiCaL*d*(fp*(1.0-nca)+jca*fcap*nca);
	ICaNa=(1.0-fICaLp)*PCaNa*PhiCaNa*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaNap*PhiCaNa*d*(fp*(1.0-nca)+jca*fcap*nca);
	ICaK=(1.0-fICaLp)*PCaK*PhiCaK*d*(f*(1.0-nca)+jca*fca*nca)+fICaLp*PCaKp*PhiCaK*d*(fp*(1.0-nca)+jca*fcap*nca);

	double xrss=1.0/(1.0+std::exp((-(v+8.337))/6.789));
	double txrf=12.98+1.0/(0.3652*std::exp((v-31.66)/3.869)+4.123e-5*std::exp((-(v-47.78))/20.38));
	double txrs=1.865+1.0/(0.06629*std::exp((v-34.70)/7.355)+1.128e-5*std::exp((-(v-29.74))/25.94));
	double Axrf=1.0/(1.0+std::exp((v+54.81)/38.21));
	double Axrs=1.0-Axrf;
	xrf=xrss-(xrss-xrf)*std::exp(-dt/txrf);
	xrs=xrss-(xrss-xrs)*std::exp(-dt/txrs);

	double xr=Axrf*xrf+Axrs*xrs;
	double rkr=1.0/(1.0+std::exp((v+55.0)/75.0))*1.0/(1.0+std::exp((v-10.0)/30.0));
	double GKr=0.046;
	if ( M_cellType == CellType::Epicardial )
	{
		GKr*=1.3;
	}
	else if ( M_cellType == CellType::MCell )
	{
		GKr*=0.8;
	}
	IKr=GKr*sqrt(ko/5.4)*xr*rkr*(v-EK);

	double xs1ss=1.0/(1.0+std::exp((-(v+11.60))/8.932));
	double txs1=817.3+1.0/(2.326e-4*std::exp((v+48.28)/17.80)+0.001292*std::exp((-(v+210.0))/230.0));
	xs1=xs1ss-(xs1ss-xs1)*std::exp(-dt/txs1);

	double xs2ss=xs1ss;
	double txs2=1.0/(0.01*std::exp((v-50.0)/20.0)+0.0193*std::exp((-(v+66.54))/31.0));
	xs2=xs2ss-(xs2ss-xs2)*std::exp(-dt/txs2);

	double KsCa=1.0+0.6/(1.0+std::pow(3.8e-5/cai,1.4));
	double GKs=0.0034;
	if ( M_cellType == CellType::Epicardial)
	{
		GKs*=1.4;
	}
	IKs=GKs*KsCa*xs1*xs2*(v-EKs);

	double xk1ss=1.0/(1.0+std::exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
	double txk1=122.2/(std::exp((-(v+127.2))/20.36)+std::exp((v+236.8)/69.33));
	xk1=xk1ss-(xk1ss-xk1)*std::exp(-dt/txk1);

	double rk1=1.0/(1.0+std::exp((v+105.8-2.6*ko)/9.493));
	double GK1=0.1908;
	if ( M_cellType == CellType::Epicardial )
	{
		GK1*=1.2;
	}
	else if ( M_cellType == CellType::MCell )
	{
		GK1*=1.3;
	}
	IK1=GK1*sqrt(ko)*rk1*xk1*(v-EK);

	double kna1=15.0;
	double kna2=5.0;
	double kna3=88.12;
	double kasymm=12.5;
	double wna=6.0e4;
	double wca=6.0e4;
	double wnaca=5.0e3;
	double kcaon=1.5e6;
	double kcaoff=5.0e3;
	double qna=0.5224;
	double qca=0.1670;
	double hca=std::exp((qca*v*F)/(R*T));
	double hna=std::exp((qna*v*F)/(R*T));
	double h1=1+nai/kna3*(1+hna);
	double h2=(nai*hna)/(kna3*h1);
	double h3=1.0/h1;
	double h4=1.0+nai/kna1*(1+nai/kna2);
	double h5=nai*nai/(h4*kna1*kna2);
	double h6=1.0/h4;
	double h7=1.0+nao/kna3*(1.0+1.0/hna);
	double h8=nao/(kna3*hna*h7);
	double h9=1.0/h7;
	double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
	double h11=nao*nao/(h10*kna1*kna2);
	double h12=1.0/h10;
	double k1=h12*cao*kcaon;
	double k2=kcaoff;
	double k3p=h9*wca;
	double k3pp=h8*wnaca;
	double k3=k3p+k3pp;
	double k4p=h3*wca/hca;
	double k4pp=h2*wnaca;
	double k4=k4p+k4pp;
	double k5=kcaoff;
	double k6=h6*cai*kcaon;
	double k7=h5*h2*wna;
	double k8=h8*h11*wna;
	double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
	double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
	double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
	double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
	double E1=x1/(x1+x2+x3+x4);
	double E2=x2/(x1+x2+x3+x4);
	double E3=x3/(x1+x2+x3+x4);
	double E4=x4/(x1+x2+x3+x4);
	double KmCaAct=150.0e-6;
	double allo=1.0/(1.0+std::pow(KmCaAct/cai,2.0));
	double zna=1.0;
	double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
	double JncxCa=E2*k2-E1*k1;
	double Gncx=0.0008;
	if ( M_cellType == CellType::Epicardial)
	{
		Gncx*=1.1;
	}
	if ( M_cellType == CellType::MCell )
	{
		Gncx*=1.4;
	}
	INaCa_i=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);

	h1=1+nass/kna3*(1+hna);
	h2=(nass*hna)/(kna3*h1);
	h3=1.0/h1;
	h4=1.0+nass/kna1*(1+nass/kna2);
	h5=nass*nass/(h4*kna1*kna2);
	h6=1.0/h4;
	h7=1.0+nao/kna3*(1.0+1.0/hna);
	h8=nao/(kna3*hna*h7);
	h9=1.0/h7;
	h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
	h11=nao*nao/(h10*kna1*kna2);
	h12=1.0/h10;
	k1=h12*cao*kcaon;
	k2=kcaoff;
	k3p=h9*wca;
	k3pp=h8*wnaca;
	k3=k3p+k3pp;
	k4p=h3*wca/hca;
	k4pp=h2*wnaca;
	k4=k4p+k4pp;
	k5=kcaoff;
	k6=h6*cass*kcaon;
	k7=h5*h2*wna;
	k8=h8*h11*wna;
	x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
	x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
	x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
	x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
	E1=x1/(x1+x2+x3+x4);
	E2=x2/(x1+x2+x3+x4);
	E3=x3/(x1+x2+x3+x4);
	E4=x4/(x1+x2+x3+x4);
	KmCaAct=150.0e-6;
	allo=1.0/(1.0+std::pow(KmCaAct/cass,2.0));
	JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
	JncxCa=E2*k2-E1*k1;
	INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);

	INaCa=INaCa_i+INaCa_ss;

	double k1p=949.5;
	double k1m=182.4;
	double k2p=687.2;
	double k2m=39.4;
	k3p=1899.0;
	double k3m=79300.0;
	k4p=639.0;
	double k4m=40.0;
	double Knai0=9.073;
	double Knao0=27.78;
	double delta=-0.1550;
	double Knai=Knai0*std::exp((delta*v*F)/(3.0*R*T));
	double Knao=Knao0*std::exp(((1.0-delta)*v*F)/(3.0*R*T));
	double Kki=0.5;
	double Kko=0.3582;
	double MgADP=0.05;
	double MgATP=9.8;
	double Kmgatp=1.698e-7;
	double H=1.0e-7;
	double eP=4.2;
	double Khp=1.698e-7;
	double Knap=224.0;
	double Kxkur=292.0;
	double P=eP/(1.0+H/Khp+nai/Knap+ki/Kxkur);
	double a1=(k1p*std::pow(nai/Knai,3.0))/(std::pow(1.0+nai/Knai,3.0)+std::pow(1.0+ki/Kki,2.0)-1.0);
	double b1=k1m*MgADP;
	double a2=k2p;
	double b2=(k2m*std::pow(nao/Knao,3.0))/(std::pow(1.0+nao/Knao,3.0)+std::pow(1.0+ko/Kko,2.0)-1.0);
	double a3=(k3p*std::pow(ko/Kko,2.0))/(std::pow(1.0+nao/Knao,3.0)+std::pow(1.0+ko/Kko,2.0)-1.0);
	double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
	double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
	double b4=(k4m*std::pow(ki/Kki,2.0))/(std::pow(1.0+nai/Knai,3.0)+std::pow(1.0+ki/Kki,2.0)-1.0);
	x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
	x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
	x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
	x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
	E1=x1/(x1+x2+x3+x4);
	E2=x2/(x1+x2+x3+x4);
	E3=x3/(x1+x2+x3+x4);
	E4=x4/(x1+x2+x3+x4);
	double zk=1.0;
	double JnakNa=3.0*(E1*a3-E2*b3);
	double JnakK=2.0*(E4*b1-E3*a1);
	double Pnak=30;
	if ( M_cellType == CellType::Epicardial )
	{
		Pnak*=0.9;
	}
	if ( M_cellType == CellType::MCell )
	{
		Pnak*=0.7;
	}
	INaK=Pnak*(zna*JnakNa+zk*JnakK);

	double xkb=1.0/(1.0+std::exp(-(v-14.48)/18.34));
	double GKb=0.003;
	if ( M_cellType == CellType::Epicardial )
	{
		GKb*=0.6;
	}
	IKb=GKb*xkb*(v-EK);

	double PNab=3.75e-10;
	INab=PNab*vffrt*(nai*std::exp(vfrt)-nao)/(std::exp(vfrt)-1.0);

	double PCab=2.5e-8;
	ICab=PCab*4.0*vffrt*(cai*std::exp(2.0*vfrt)-0.341*cao)/(std::exp(2.0*vfrt)-1.0);

	double GpCa=0.0005;
	IpCa=GpCa*cai/(0.0005+cai);

	//variables[1] = nai;
	//variables[2] = nass;
	//variables[3] = ki;
	//variables[4] = kss;
	//variables[5] = cai;
	//variables[6] = cass;
	//variables[7] = cansr;
	//variables[8] = cajsr;
	variables[9] = m;
	variables[10] = hf;
	variables[11]= hs;
	variables[12]= j;
	variables[13]= hsp;
	variables[14]= jp;
	variables[15]= mL;
	variables[16]= hL;
	variables[17]= hLp;
	variables[18]= a;
	variables[19]= iF;
	variables[20]= iS;
	variables[21]= ap;
	variables[22]= iFp;
	variables[23]= iSp;
	variables[24]= d;
	variables[25]= ff;
	variables[26]= fs;
	variables[27]= fcaf;
	variables[28]= fcas;
	variables[29]= jca;
	variables[30]= nca;
	variables[31]= ffp;
	variables[32]= fcafp;
	variables[33]= xrf;
	variables[34]= xrs;
	variables[35]= xs1;
	variables[36]= xs2;
	variables[37]= xk1;
	//variables[38]= Jrelnp;
	//variables[39]= Jrelp;
	//variables[40]= CaMKt;
}


void
ORd::FBC(std::vector<double>& variables, double dt)
{
    double nai   = variables[1];
    double nass  = variables[2];
    double ki    = variables[3];
    double kss   = variables[4];
    double cai   = variables[5];
    double cass  = variables[6];
    double cansr = variables[7];
    double cajsr = variables[8];

    double Jrelnp= variables[38];
    double Jrelp = variables[39];
    double CaMKt = variables[40];

    double CaMKb=CaMKo*(1.0-CaMKt)/(1.0+KmCaM/cass);
    CaMKa=CaMKb+CaMKt;
    CaMKt+=dt*(aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt);

    JdiffNa=(nass-nai)/2.0;
    JdiffK=(kss-ki)/2.0;
    Jdiff=(cass-cai)/0.2;

    double bt=4.75;
    double a_rel=0.5*bt;
    double Jrel_inf=a_rel*(-ICaL)/(1.0+std::pow(1.5/cajsr,8.0));
    if (M_cellType == CellType::MCell )
    {
        Jrel_inf*=1.7;
    }
    double tau_rel=bt/(1.0+0.0123/cajsr);
    if (tau_rel<0.005)
    {
        tau_rel=0.005;
    }
    Jrelnp=Jrel_inf-(Jrel_inf-Jrelnp)*std::exp(-dt/tau_rel);

    double btp=1.25*bt;
    double a_relp=0.5*btp;
    double Jrel_infp=a_relp*(-ICaL)/(1.0+std::pow(1.5/cajsr,8.0));
    if (M_cellType==CellType::MCell)
    {
        Jrel_infp*=1.7;
    }
    double tau_relp=btp/(1.0+0.0123/cajsr);
    if (tau_relp<0.005)
    {
        tau_relp=0.005;
    }
    Jrelp=Jrel_infp-(Jrel_infp-Jrelp)*std::exp(-dt/tau_relp);

    double fJrelp=(1.0/(1.0+KmCaMK/CaMKa));
    Jrel=(1.0-fJrelp)*Jrelnp+fJrelp*Jrelp;

    double Jupnp=0.004375*cai/(cai+0.00092);
    double Jupp=2.75*0.004375*cai/(cai+0.00092-0.00017);
    if (M_cellType == CellType::Epicardial )
    {
        Jupnp*=1.3;
        Jupp*=1.3;
    }
    double fJupp=(1.0/(1.0+KmCaMK/CaMKa));
    Jleak=0.0039375*cansr/15.0;
    Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;

    Jtr=(cansr-cajsr)/100.0;

    nai+=dt*(-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo);
    nass+=dt*(-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa);

    ki+=dt*(-(Ito+IKr+IKs+IK1+IKb+Ist-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo);
    kss+=dt*(-(ICaK)*Acap/(F*vss)-JdiffK);

    double Bcai;
    if (M_cellType == CellType::Epicardial)
    {
        Bcai=1.0/(1.0+1.3*cmdnmax*kmcmdn/std::pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/std::pow(kmtrpn+cai,2.0));
    }
    else
    {
        Bcai=1.0/(1.0+cmdnmax*kmcmdn/std::pow(kmcmdn+cai,2.0)+trpnmax*kmtrpn/std::pow(kmtrpn+cai,2.0));
    }
    cai+=dt*(Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo));

    double Bcass=1.0/(1.0+BSRmax*KmBSR/std::pow(KmBSR+cass,2.0)+BSLmax*KmBSL/std::pow(KmBSL+cass,2.0));
    cass+=dt*(Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff));

    cansr+=dt*(Jup-Jtr*vjsr/vnsr);

    double Bcajsr=1.0/(1.0+csqnmax*kmcsqn/std::pow(kmcsqn+cajsr,2.0));
    cajsr+=dt*(Bcajsr*(Jtr-Jrel));

    variables[1] = nai;
    variables[2] = nass;
    variables[3] = ki;
    variables[4] = kss;
    variables[5] = cai;
    variables[6] = cass;
    variables[7] = cansr;
    variables[8] = cajsr;


    variables[38]= Jrelnp;
    variables[39]= Jrelp;
    variables[40]= CaMKt;
}

} /* namespace BeatIt */
