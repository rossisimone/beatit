/*
 ============================================================================

 .______    _______     ___   .___________.    __  .___________.
 |   _  \  |   ____|   /   \  |           |   |  | |           |
 |  |_)  | |  |__     /  ^  \ `---|  |----`   |  | `---|  |----`
 |   _  <  |   __|   /  /_\  \    |  |        |  |     |  |     
 |  |_)  | |  |____ /  _____  \   |  |        |  |     |  |     
 |______/  |_______/__/     \__\  |__|        |__|     |__|     
 
 BeatIt - code for cardiovascular simulations
 Copyrigdt (C) 2016 Simone Rossi

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <dttp://www.gnu.org/licenses/>.
 ============================================================================
 */

/**
 * \file TP06.cpp
 *
 * \author srossi
 *
 * \version 0.0
 *
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Sep 2, 2016
 *
 */

#include <cmath>
#include <fstream>
#include "Electrophysiology/IonicModels/TP06.hpp"

namespace BeatIt
{


IonicModel* createTP06()
{
	return new TP06();
}

void
TP06::setCellType(CellType type)
{
    super::setCellType(type);
    selectParameters(M_cellType);
}

void
TP06::selectParameters(CellType type)
{
    switch(M_cellType)
    {
        case CellType::Epicardial:
        {
            Gks=0.392;
            Gto=0.294;
            break;
        }
        case CellType::Endocardial:
        {
            Gks=0.392;
            Gto=0.073;
            break;
        }
        case CellType::MCell:
        {
            Gks=0.098;
            Gto=0.294;
            break;
        }
    }
}

TP06::TP06()
  : super(20,0, "TP06", CellType::MCell)
{
    // Without potential

    //Cai=0.00007;
    M_variablesNames[0] = "Cai";
    //CaSR_init=1.3
    M_variablesNames[1] = "CaSR";
    //CaSS_init=0.00007
    M_variablesNames[2] = "CaSS";
    //Nai_init=7.67
    M_variablesNames[3] = "Nai";
    //Ki_init=138.3
    M_variablesNames[4] = "Ki";
    //M= 0.
    M_variablesNames[5] = "M";
    //H= 0.75
    M_variablesNames[6] = "H";
    //J= 0.75
    M_variablesNames[7] = "J";
    //Xr1= 0.
    M_variablesNames[8] = "Xr1";
    //Xr2= 1.
    M_variablesNames[9] = "Xr2";
    //Xs= 0.
    M_variablesNames[10] = "Xs";
    // S= 1.
    M_variablesNames[11] = "S";
    //R= 0.
    M_variablesNames[12] = "R";
    //D= 0.
    M_variablesNames[13] = "D";
    //F= 1.
    M_variablesNames[14] = "F";
    //F2=1.
    M_variablesNames[15] = "F2";
    //FCass= 1.
    M_variablesNames[16] = "FCass";
    //RR= 1.
    M_variablesNames[17] = "RR";
    //OO= 0.
    M_variablesNames[18] = "OO";

    selectParameters(M_cellType);
}

void
TP06::initialize(std::vector<double>& variables)
{
    // V
	variables[0] = -86.2;
    //Cai=0.00007
    variables[1] = 0.00007;
    //CaSR_init=1.3
    variables[2] = 1.3;
    //CaSS_init=0.00007
    variables[3] = 0.00007;
    //Nai_init=7.67
    variables[4] = 7.67;
    //Ki_init=138.3
    variables[5] = 138.3;
    //M= 0.
    variables[6] = 0.0;
    //H= 0.75
    variables[7] = 0.75;
    //J= 0.75
    variables[8] = 0.75;
    //Xr1= 0.
    variables[9] = 0.0;
    //Xr2= 1.
    variables[10] = 1.0;
    //Xs= 0.
    variables[11] = 0.0;
    //S= 1.
    variables[12] = 1.0;
    //R= 0.
    variables[13] = 0.0;
    //D= 0.
    variables[14] = 0.0;
    //F= 1.
    variables[15] = 1.0;
    //F2=1.
    variables[16] = 1.0;
    //FCass= 1.
    variables[17] = 1.0;
    //RR= 1.
    variables[18] = 1.0;
    //OO= 0.
    variables[19] = 0.0;
}


//! Update all the variables in the ionic model
/*!
 *  \param [in] variables Vector containing the local value of all variables
 *  \param [in] dt        Timestep
 */
void
TP06::updateVariables(std::vector<double>& variables, double appliedCurrent, double dt)
{
	// For compatibility  with the original code where the applied stimulus in opposite
	Istim = appliedCurrent;
    step(variables, dt);
}

//! Evaluate total ionic current for the computation of the potential
/*!
 *  \param [in] variables Vector containing the local value of all variables
 *  \param [in] appliedCurrent value of the applied current
 *  \param [in] dt        Timestep
 */
double
TP06::evaluateIonicCurrent(std::vector<double>& variables, double appliedCurrent, double dt)
{
	return Itot;

}

double
TP06::evaluateIonicCurrentTimeDerivative( std::vector<double>& variables,
                                     std::vector<double>& old_variables,
                                     double dt,
                                     double h )
{
    double svolt = variables[0];
    double Q = old_variables[0];
    double dIdV = dItot;

    double Cai   = variables[1];
    double CaSR  = variables[2];
    double CaSS  = variables[3];
    double Nai   = variables[4];
    double Ki    = variables[5];
    double sm    = variables[6];
    double sh    = variables[7];
    double sj    = variables[8];
    double sxr1  = variables[9];
    double sxr2  = variables[10];
    double sxs   = variables[11];
    double ss    = variables[12];
    double sr    = variables[13];
    double sd    = variables[14];
    double sf    = variables[15];
    double sf2   = variables[16];
    double sfcass= variables[17];
    double dCai   = (variables[1]- old_variables[1])/dt;
    double dCaSR  = (variables[2]- old_variables[2])/dt;
    double dCaSS  = (variables[3]- old_variables[3])/dt;
    double dNai   = (variables[4]- old_variables[4])/dt;
    double dKi    = (variables[5]- old_variables[5])/dt;
    double dsm    = (variables[6]- old_variables[6])/dt;
    double dsh    = (variables[7]- old_variables[7])/dt;
    double dsj    = (variables[8]- old_variables[8])/dt;
    double dsxr1  = (variables[9]- old_variables[9])/dt;
    double dsxr2  = (variables[10]- old_variables[10])/dt;
    double dsxs   = (variables[11]- old_variables[11])/dt;
    double dss    = (variables[12]- old_variables[12])/dt;
    double dsr    = (variables[13]- old_variables[13])/dt;
    double dsd    = (variables[14]- old_variables[14])/dt;
    double dsf    = (variables[15]- old_variables[15])/dt;
    double dsf2   = (variables[16]- old_variables[16])/dt;
    double dsfcass= (variables[17]- old_variables[17])/dt;


    //Compute currents
    double dINa=3*GNa*sm*sm*dsm*sh*sj*(svolt-Ena)
               +GNa*sm*sm*sm*dsh*sj*(svolt-Ena)
               +GNa*sm*sm*sm*sh*dsj*(svolt-Ena);
    double aux1_dICaL = (std::exp(2*(svolt-15)*F/(R*T))-1.); // denominator
    double dICaL=GCaL*sd*sf*sf2*sfcass*4*(svolt-15)*(F*F/(R*T))*
      (0.25*std::exp(2*(svolt-15)*F/(R*T))*dCaSS)/ aux1_dICaL // dCaSS
      + GCaL*dsd*sf*sf2*sfcass*4*(svolt-15)*(F*F/(R*T))*
      (0.25*std::exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/ aux1_dICaL// dsd
    + GCaL*sd*dsf*sf2*sfcass*4*(svolt-15)*(F*F/(R*T))*
    (0.25*std::exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/ aux1_dICaL// dsf
    + GCaL*sd*sf*dsf2*sfcass*4*(svolt-15)*(F*F/(R*T))*
    (0.25*std::exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/ aux1_dICaL//dsf2
      + GCaL*sd*sf*sf2*dsfcass*4*(svolt-15)*(F*F/(R*T))*
      (0.25*std::exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/ aux1_dICaL; //dsfcass
    double dIto=Gto*dsr*ss*(svolt-Ek)
               +Gto*sr*dss*(svolt-Ek);
    double dIKr=Gkr*sqrt(Ko/5.4)*dsxr1*sxr2*(svolt-Ek)
               +Gkr*sqrt(Ko/5.4)*sxr1*dsxr2*(svolt-Ek);
    double dIKs=2*Gks*sxs*dsxs*(svolt-Eks);

    double dINaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
      (1./(1+ksat*std::exp((n-1)*svolt*F/(R*T))))*
      (std::exp(n*svolt*F/(R*T))*3*Nai*Nai*dNai*Cao-
       std::exp((n-1)*svolt*F/(R*T))*Nao*Nao*Nao*dCai*2.5);
    //INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    double aux_dINak = dNai/(Nai+KmNa) - Nai/(Nai+KmNa) /(Nai+KmNa) * dNai;
    double dINaK=knak*(Ko/(Ko+KmK))*aux_dINak*rec_iNaK;
    //IpCa=GpCa*Cai/(KpCa+Cai);
    double dIpCa=GpCa*(dCai/(KpCa+Cai)-Cai/(KpCa+Cai)/(KpCa+Cai)*dCai);
    //IpK=GpK*rec_ipK*(svolt-Ek);
    //IbNa=GbNa*(svolt-Ena);
    //IbCa=GbCa*(svolt-Eca);
    double dI = dIdV * Q + dINa + dICaL + dIto + dIKr + dIKs + dINaCa + dINaK + dIpCa;
    return dI;

}

void
TP06::initializeSaveData(std::ostream& output)
{
	// time -  0
	output << "time v ";
	//  1 - 10
	output << "Cai CaSR CaSS Nai Ki M H J Xr1 Xr2 ";
	// 11 - 19
	output << "Xs S R D F F2 FCass RR OO\n";
}

void
TP06::step(std::vector<double>& variables, double dt)
{
	double& svolt = variables[0];
	double& Cai   = variables[1];
	double& CaSR  = variables[2];
	double& CaSS  = variables[3];
	double& Nai   = variables[4];
	double& Ki    = variables[5];
	double& sm    = variables[6];
	double& sh    = variables[7];
	double& sj    = variables[8];
	double& sxr1  = variables[9];
	double& sxr2  = variables[10];
	double& sxs   = variables[11];
	double& ss    = variables[12];
	double& sr    = variables[13];
	double& sd    = variables[14];
	double& sf    = variables[15];
	double& sf2   = variables[16];
	double& sfcass= variables[17];
	double& sRR   = variables[18];
	double& sOO   = variables[19];
	double& svolt2= Volt2;
	double& sItot = Itot;

    //Needed to compute currents
    Ek=RTONF*(std::log((Ko/Ki)));
    Ena=RTONF*(std::log((Nao/Nai)));
    Eks=RTONF*(std::log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    Eca=0.5*RTONF*(std::log((Cao/Cai)));
    Ak1=0.1/(1.+std::exp(0.06*(svolt-Ek-200)));
    Bk1=(3.*std::exp(0.0002*(svolt-Ek+100))+
     std::exp(0.1*(svolt-Ek-10)))/(1.+std::exp(-0.5*(svolt-Ek)));
    rec_iK1=Ak1/(Ak1+Bk1);
    rec_iNaK=(1./(1.+0.1245*std::exp(-0.1*svolt*F/(R*T))+0.0353*std::exp(-svolt*F/(R*T))));
    rec_ipK=1./(1.+std::exp((25-svolt)/5.98));


    //Compute currents
    INa=GNa*sm*sm*sm*sh*sj*(svolt-Ena);
    ICaL=GCaL*sd*sf*sf2*sfcass*4*(svolt-15)*(F*F/(R*T))*
      (0.25*std::exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/(std::exp(2*(svolt-15)*F/(R*T))-1.);
    Ito=Gto*sr*ss*(svolt-Ek);
    IKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(svolt-Ek);
    IKs=Gks*sxs*sxs*(svolt-Eks);
    IK1=GK1*rec_iK1*(svolt-Ek);
    INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
      (1./(1+ksat*std::exp((n-1)*svolt*F/(R*T))))*
      (std::exp(n*svolt*F/(R*T))*Nai*Nai*Nai*Cao-
       std::exp((n-1)*svolt*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
    INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    IpCa=GpCa*Cai/(KpCa+Cai);
    IpK=GpK*rec_ipK*(svolt-Ek);
    IbNa=GbNa*(svolt-Ena);
    IbCa=GbCa*(svolt-Eca);


    //Determine total current
    Itot = IKr   +
           IKs   +
           IK1   +
           Ito   +
           INa   +
           IbNa  +
           ICaL  +
           IbCa  +
           INaK  +
           INaCa +
           IpCa  +
           IpK;/*   +
           Istim;*/

    //Compute currents derivatives
    double dINa=GNa*sm*sm*sm*sh*sj;
    double dICaL=(2.0*CaSS*F*F*F*GCaL*sd*sf*sf2*sfcass*std::exp((F*(2.0*svolt - 30.0))/(R*T))
                * (svolt - 15.0))/(R*R*T*T*(std::exp((F*(2.0*svolt - 30.0))/(R*T)) - 1.0))
                - (4.0*F*F*GCaL*sd*sf*sf2*sfcass
                * (Cao - 0.25*CaSS*std::exp((F*(2.0*svolt - 30.0))/(R*T))))
                / (R*T*(std::exp((F*(2.0*svolt - 30.0))/(R*T)) - 1.0))
                + (8.0*F*F*F*GCaL*sd*sf*sf2*sfcass*std::exp((F*(2.0*svolt - 30.0))/(R*T))
                * (svolt - 15.0) * (Cao - 0.25*CaSS*std::exp((F*(2.0*svolt - 30.0))/(R*T))))
                /(R*R*T*T*(std::exp((F*(2.0*svolt - 30.0))/(R*T)) - 1.0)*(std::exp((F*(2.0*svolt - 30.0))/(R*T)) - 1.0));
    double dIto=Gto*sr*ss;
    double dIKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2;
    double dIKs=Gks*sxs*sxs;

    double dAk1= -(0.006*std::exp(0.06*svolt - 0.06*Ek - 12.0))
               / (std::exp(0.06*svolt - 0.06*Ek - 12.0) + 1.0)
               / (std::exp(0.06*svolt - 0.06*Ek - 12.0) + 1.0);
    double dBk1= (0.0006*std::exp(0.0002*svolt - 0.0002*Ek + 0.02)
               + 0.1*std::exp(0.1*svolt - 0.1*Ek - 1.0))
               / ( std::exp(0.5*Ek - 0.5*svolt) + 1.0)
               + (0.5*std::exp(0.5*Ek - 0.5*svolt)
               * (3.0*std::exp(0.0002*svolt - 0.0002*Ek + 0.02)
               + std::exp(0.1*svolt - 0.1*Ek - 1.0)))
               / (std::exp(0.5*Ek - 0.5*svolt) + 1.0)
               / (std::exp(0.5*Ek - 0.5*svolt) + 1.0);
    double drec_iK1=dAk1/(Ak1+Bk1) - Ak1/(Ak1+Bk1)/(Ak1+Bk1)*(dAk1+dBk1);
    double drec_iNaK= ((0.01245*F*std::exp(-(0.1*F*svolt)/(R*T)))/(R*T)
                    + (0.0353*F*std::exp(-(1.0*F*svolt)/(R*T)))/(R*T))
                    / (0.1245*std::exp(-(0.1*F*svolt)/(R*T))+0.0353*exp(-(1.0*F*svolt)/(R*T))+1.0)
                    / (0.1245*std::exp(-(0.1*F*svolt)/(R*T))+0.0353*exp(-(1.0*F*svolt)/(R*T))+1.0);
    double drec_ipK=  (std::exp((25-svolt)/5.98)) / 5.98
                   / (1.+std::exp((25-svolt)/5.98))
                   /(1.+std::exp((25-svolt)/5.98));

    double dIK1=GK1*drec_iK1*(svolt-Ek)+GK1*rec_iK1;
    double dINaCa=(1.0*knaca*((1.0*Cao*F*n*std::exp((F*n*svolt)/(R*T))*Nai*Nai*Nai)/(R*T)
                 - (2.5*Cai*F*Nao*Nao*Nao*std::exp((F*svolt*(n - 1.0))/(R*T))*(n - 1.0))/(R*T)))
                 / ((KmNai*KmNai*KmNai + Nao*Nao*Nao)
                 * (ksat*std::exp((F*svolt*(n - 1.0))/(R*T)) + 1.0)*(Cao + KmCa))
                 - (F*knaca*ksat*std::exp((F*svolt*(n - 1.0))/(R*T))*(n - 1.0)
                 * (1.0*Cao*std::exp((F*n*svolt)/(R*T))*Nai*Nai*Nai
                 - 2.5*Cai*std::exp((F*svolt*(n - 1.0))/(R*T))*Nao*Nao*Nao))
                 / (R*T*(KmNai*KmNai*KmNai + Nao*Nao*Nao)
                 * (ksat*std::exp((F*svolt*(n - 1.0))/(R*T)) + 1.0)
                 * (ksat*std::exp((F*svolt*(n - 1.0))/(R*T)) + 1.0)*(Cao + KmCa));
    double dINaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*drec_iNaK;
    double dIpCa=0.0;
    double dIpK=GpK*drec_ipK*(svolt-Ek)+GpK*rec_ipK;
    double dIbNa=GbNa;
    double dIbCa=GbCa;


    //Determine total current
    dItot = dIKr   +
            dIKs   +
            dIK1   +
            dIto   +
            dINa   +
            dIbNa  +
            dICaL  +
            dIbCa  +
            dINaK  +
            dINaCa +
            dIpCa  +
            dIpK;
    //update concentrations
    kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR)));
    k1=k1_/kCaSR;
    k2=k2_*kCaSR;
    dRR=k4*(1-sRR)-k2*CaSS*sRR;
    sRR+=dt*dRR;
    sOO=k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);


    Irel=Vrel*sOO*(CaSR-CaSS);
    Ileak=Vleak*(CaSR-Cai);
    Iup=Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));
    Ixfer=Vxfer*(CaSS-Cai);


    CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
    dCaSR=dt*(Iup-Irel-Ileak);
    bjsr=Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
    cjsr=Kbufsr*(CaCSQN+dCaSR+CaSR);
    CaSR=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;


    CaSSBuf=Bufss*CaSS/(CaSS+Kbufss);
    dCaSS=dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE));
    bcss=Bufss-CaSSBuf-dCaSS-CaSS+Kbufss;
    ccss=Kbufss*(CaSSBuf+dCaSS+CaSS);
    CaSS=(sqrt(bcss*bcss+4*ccss)-bcss)/2;


    CaBuf=Bufc*Cai/(Cai+Kbufc);
    dCai=dt*((-(IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer);
    bc=Bufc-CaBuf-dCai-Cai+Kbufc;
    cc=Kbufc*(CaBuf+dCai+Cai);
    Cai=(sqrt(bc*bc+4*cc)-bc)/2;


    dNai=-(INa+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE;
    Nai+=dt*dNai;

    dKi=-(Istim+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;
    Ki+=dt*dKi;



    //compute steady state values and time constants
    AM=1./(1.+std::exp((-60.-svolt)/5.));
    BM=0.1/(1.+std::exp((svolt+35.)/5.))+0.10/(1.+std::exp((svolt-50.)/200.));
    TAU_M=AM*BM;
    M_INF=1./((1.+std::exp((-56.86-svolt)/9.03))*(1.+std::exp((-56.86-svolt)/9.03)));
    if (svolt>=-40.)
    {
        AH_1=0.;
        BH_1=(0.77/(0.13*(1.+std::exp(-(svolt+10.66)/11.1))));
        TAU_H= 1.0/(AH_1+BH_1);
    }
    else
    {
        AH_2=(0.057*std::exp(-(svolt+80.)/6.8));
        BH_2=(2.7*std::exp(0.079*svolt)+(3.1e5)*std::exp(0.3485*svolt));
        TAU_H=1.0/(AH_2+BH_2);
    }
    H_INF=1./((1.+std::exp((svolt+71.55)/7.43))*(1.+std::exp((svolt+71.55)/7.43)));
    if(svolt>=-40.)
    {
        AJ_1=0.;
        BJ_1=(0.6*std::exp((0.057)*svolt)/(1.+std::exp(-0.1*(svolt+32.))));
        TAU_J= 1.0/(AJ_1+BJ_1);
    }
    else
    {
        AJ_2=(((-2.5428e4)*std::exp(0.2444*svolt)-(6.948e-6)*
          std::exp(-0.04391*svolt))*(svolt+37.78)/
             (1.+std::exp(0.311*(svolt+79.23))));
        BJ_2=(0.02424*std::exp(-0.01052*svolt)/(1.+std::exp(-0.1378*(svolt+40.14))));
        TAU_J= 1.0/(AJ_2+BJ_2);
    }
    J_INF=H_INF;

    Xr1_INF=1./(1.+std::exp((-26.-svolt)/7.));
    axr1=450./(1.+std::exp((-45.-svolt)/10.));
    bxr1=6./(1.+std::exp((svolt-(-30.))/11.5));
    TAU_Xr1=axr1*bxr1;
    Xr2_INF=1./(1.+std::exp((svolt-(-88.))/24.));
    axr2=3./(1.+std::exp((-60.-svolt)/20.));
    bxr2=1.12/(1.+std::exp((svolt-60.)/20.));
    TAU_Xr2=axr2*bxr2;

    Xs_INF=1./(1.+std::exp((-5.-svolt)/14.));
    Axs=(1400./(sqrt(1.+std::exp((5.-svolt)/6))));
    Bxs=(1./(1.+std::exp((svolt-35.)/15.)));
    TAU_Xs=Axs*Bxs+80;


    switch(M_cellType)
    {
      case CellType::Endocardial:
      {
          R_INF=1./(1.+std::exp((20-svolt)/6.));
          S_INF=1./(1.+std::exp((svolt+28)/5.));
          TAU_R=9.5*std::exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
          TAU_S=1000.*std::exp(-(svolt+67)*(svolt+67)/1000.)+8.;
          break;
      }
      case CellType::Epicardial:
      {
          R_INF=1./(1.+std::exp((20-svolt)/6.));
          S_INF=1./(1.+std::exp((svolt+20)/5.));
          TAU_R=9.5*std::exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
          TAU_S=85.*std::exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+std::exp((svolt-20.)/5.))+3.;
          break;
      }
      default:
      case CellType::MCell:
      {
          R_INF=1./(1.+std::exp((20-svolt)/6.));
          S_INF=1./(1.+std::exp((svolt+20)/5.));
          TAU_R=9.5*std::exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
          TAU_S=85.*std::exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+std::exp((svolt-20.)/5.))+3.;
          break;
      }
    }

    D_INF=1./(1.+std::exp((-8-svolt)/7.5));
    Ad=1.4/(1.+std::exp((-35-svolt)/13))+0.25;
    Bd=1.4/(1.+std::exp((svolt+5)/5));
    Cd=1./(1.+std::exp((50-svolt)/20));
    TAU_D=Ad*Bd+Cd;
    F_INF=1./(1.+std::exp((svolt+20)/7));
    Af=1102.5*std::exp(-(svolt+27)*(svolt+27)/225);
    Bf=200./(1+std::exp((13-svolt)/10.));
    Cf=(180./(1+std::exp((svolt+30)/10)))+20;
    TAU_F=Af+Bf+Cf;
    F2_INF=0.67/(1.+std::exp((svolt+35)/7))+0.33;
    Af2=600*std::exp(-(svolt+25)*(svolt+25)/170);
    Bf2=31/(1.+std::exp((25-svolt)/10));
    Cf2=16/(1.+std::exp((svolt+30)/10));
    TAU_F2=Af2+Bf2+Cf2;
    FCaSS_INF=0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;
    TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;

    sm = M_INF-(M_INF-sm)*std::exp(-dt/TAU_M);
    sh = H_INF-(H_INF-sh)*std::exp(-dt/TAU_H);
    sj = J_INF-(J_INF-sj)*std::exp(-dt/TAU_J);
    sxr1 = Xr1_INF-(Xr1_INF-sxr1)*std::exp(-dt/TAU_Xr1);
    sxr2 = Xr2_INF-(Xr2_INF-sxr2)*std::exp(-dt/TAU_Xr2);
    sxs = Xs_INF-(Xs_INF-sxs)*std::exp(-dt/TAU_Xs);
    ss= S_INF-(S_INF-ss)*std::exp(-dt/TAU_S);
    sr= R_INF-(R_INF-sr)*std::exp(-dt/TAU_R);
    sd = D_INF-(D_INF-sd)*std::exp(-dt/TAU_D);
    sf =F_INF-(F_INF-sf)*std::exp(-dt/TAU_F);
    sf2 =F2_INF-(F2_INF-sf2)*std::exp(-dt/TAU_F2);
    sfcass =FCaSS_INF-(FCaSS_INF-sfcass)*std::exp(-dt/TAU_FCaSS);

}

} /* namespace BeatIt */
