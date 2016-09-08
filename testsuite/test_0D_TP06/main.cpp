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
 * \file main.cpp
 *
 * \brief This example uses the original TP source code
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




/* Header file, forms a generalized header file containing headers of libraries that are used,
   function definitions and the cell type specification */

#include <stdio.h>
#include <cmath>
//#include <stdlib>
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;
/*------------------------------------------------------------------------------
FLAG TO CHOOSE BETWEEN EPICARDIAL ENDOCARDIAL AND MIDMYOCARDIAL CELL TYPES
------------------------------------------------------------------------------*/
#define EPI
//#define ENDO
//#define MCELL



class Variables
{
public:
  //voltage + backup
  double   Volt;
  double   Volt2;
  double   Cai;
  double   CaSR;
  double   CaSS;
  double   Nai;
  double   Ki;


  //states of voltage and time dependent gates
  //INa
  double   M;
  double   H;
  double   J;
  //IKr
  //IKr1
  double   Xr1;
  //IKr2
  double   Xr2;
  //IKs
  double   Xs;
  //Ito1
  double   R;
  double   S;
  //ICa
  double   D;
  double   F;
  double   F2;
  double   FCass;
  //Irel
  double   RR;
  double   OO;

  //total current
  double   Itot;


public:
  Variables(double V_init,double Cai_init,double CaSR_init,double CaSS_init,double Nai_init,double Ki_init);
  void writebackup(double *t,char *despath);
};


Variables::Variables(double V_init,double Cai_init,double CaSR_init,double CaSS_init,double Nai_init,double Ki_init)
{
  Volt=V_init;
  Volt2=V_init;
  Cai=Cai_init;
  CaSR=CaSR_init;
  CaSS=CaSS_init;
  Nai=Nai_init;
  Ki=Ki_init;
  M= 0.;
  H= 0.75;
  J= 0.75;
  Xr1= 0.;
  Xr2= 1.;
  Xs= 0.;
  R= 0.;
  S= 1.;
  D= 0.;
  F= 1.;
  F2=1.;
  FCass= 1.;
  RR= 1.;
  OO= 0.;

  printf("Variables initialized\n");
}



void Variables::writebackup(double *t,char *despath)
{
  static char filename[300];

  sprintf(filename,"%s%s",despath,"/PointBackupData");

  ofstream oo(filename,ios::app);
  if(!oo)
    {
      printf("cannot open file %s\n",filename);
      exit(1);
    }

  oo << *t << "\t";
  oo << Volt<< "\t";
  oo << Volt2 << "\t";
  oo << Cai << "\t";
  oo << CaSR << "\t";
  oo << CaSS << "\t";
  oo << Nai << "\t";
  oo << Ki << "\t";
  oo << M << "\t";
  oo << H << "\t";
  oo << J << "\t";
  oo << Xr1 << "\t";
  oo << Xr2 << "\t";
  oo << Xs << "\t";
  oo << S << "\t";
  oo << R << "\t";
  oo << D << "\t";
  oo << F << "\t";
  oo << F2 << "\t";
  oo << FCass << "\t";
  oo << RR << "\t";
  oo << OO << "\t";
  oo << Itot;
  oo << endl;
  oo.close();

}



void Step(Variables *V,double HT,char *despath,double *t,int step,
          double Istim);





/*-----------------------------------------------------------------------------
  FLAG TO CHOOSE BETWEEN DYNAMIC AND S1S2 RESTITUTION PROTOCOL
-----------------------------------------------------------------------------*/
#define DYNRESTPROTOCOL
//#define  S1S2RESTPROTOCOL

/*-----------------------------------------------------------------------------
  ELECTROPHYSIOLOGICAL PARAMETERS:
-----------------------------------------------------------------------------*/

//External concentrations
double Ko=5.4;
double Cao=2.0;
double Nao=140.0;

//Intracellular volumes
double Vc=0.016404;
double Vsr=0.001094;
double Vss=0.00005468;

//Calcium buffering dynamics
double Bufc=0.2;
double Kbufc=0.001;
double Bufsr=10.;
double Kbufsr=0.3;
double Bufss=0.4;
double Kbufss=0.00025;

//Intracellular calcium flux dynamics
double Vmaxup=0.006375;
double Kup=0.00025;
double Vrel=0.102;//40.8;
double k1_=0.15;
double k2_=0.045;
double k3=0.060;
double k4=0.005;//0.000015;
double EC=1.5;
double maxsr=2.5;
double minsr=1.;
double Vleak=0.00036;
double Vxfer=0.0038;



//Constants
double R=8314.472;
double F=96485.3415;
double T=310.0;
double RTONF=(R*T)/F;

//Cellular capacitance
double CAPACITANCE=0.185;

//Parameters for currents
//Parameters for IKr
double Gkr=0.153;
//Parameters for Iks
double pKNa=0.03;
#ifdef EPI
double Gks=0.392;
#endif
#ifdef ENDO
double Gks=0.392;
#endif
#ifdef MCELL
double Gks=0.098;
#endif
//Parameters for Ik1
double GK1=5.405;
//Parameters for Ito
#ifdef EPI
double Gto=0.294;
#endif
#ifdef ENDO
double Gto=0.073;
#endif
#ifdef MCELL
double Gto=0.294;
#endif
//Parameters for INa
double GNa=14.838;
//Parameters for IbNa
double GbNa=0.00029;
//Parameters for INaK
double KmK=1.0;
double KmNa=40.0;
double knak=2.724;
//Parameters for ICaL
double GCaL=0.00003980;
//Parameters for IbCa
double GbCa=0.000592;
//Parameters for INaCa
double knaca=1000;
double KmNai=87.5;
double KmCa=1.38;
double ksat=0.1;
double n=0.35;
//Parameters for IpCa
double GpCa=0.1238;
double KpCa=0.0005;
//Parameters for IpK;
double GpK=0.0146;


/*------------------------------------------------------------------------------
                PARAMETER FOR INTEGRATION
------------------------------------------------------------------------------*/
//timestep (ms)
double  HT =0.02;

/*-----------------------------------------------------------------------------
                PARAMETERS FOR INITIAL CONDITIONS
------------------------------------------------------------------------------*/
//Initial values of state variables
double V_init=-86.2;
double Cai_init=0.00007;
double CaSR_init=1.3;
double CaSS_init=0.00007;
double Nai_init=7.67;
double Ki_init=138.3;

/*--------------------------------------- ------------------------------------
             PARAMETER FOR SIMULATION DURATION
  ---------------------------------------------------------------------------*/
//duration of the simulation
double STOPTIME=100000;

/*-----------------------------------------------------------------------------
  PARAMETERS FOR STIMULATION PROTOCOLS
-----------------------------------------------------------------------------*/

#ifdef DYNRESTPROTOCOL
int i_low=0,i_high=1;
int j_low=0,j_high=1;
double stimduration=1.0;
double stimstrength=-38;
double period=1000;
double sum=period*100.;
double tbegin=0;
double tend=tbegin+stimduration;
#endif


#ifdef S1S2RESTPROTOCOL
int i_low=0,i_high=1;
int j_low=0,j_high=1;
double stimduration=1.;
double stimstrength=-52;
double tbegin=0;
double tend=tbegin+stimduration;
int counter=1;
double dia=5000;
double basicperiod=1000.;
double basicapd=274;
int repeats=10;
#endif

/*----------------------------------------------------------------------------
                            OTHER PARAMETERS
  ---------------------------------------------------------------------------*/

//destination path to put in output files
char despath[200];

/*---------------------------------------------------------------------------*/


#define FALSE 0
#define TRUE 1



void Start(int argc, char *argv[])
{
  if(argc<2)
    {
      printf("type: model destinationpath\n");
      exit(1);
    }
  else
    {
      std::cout << "copying string ... " << std::flush;
      strcpy(despath,argv[1]);
      std::cout << "done." << std::endl;
    }
}

int main(int argc, char *argv[])
{
  static double time=0;
  int step;
  double Istim;

  Start(argc,argv);

  Variables Var(V_init,Cai_init,CaSR_init,CaSS_init,Nai_init,Ki_init);


  for(step=0;time<=STOPTIME;step++)
    {
      time+=HT;

#ifdef DYNRESTPROTOCOL
      if(time>sum)
    {
      if(period>4000)
        {
          period=period-1000;
          sum=sum+period*30;
        }
      else if(period>3000)
        {
          period=period-1000;
          sum=sum+period*30;
        }
      else if(period>2000)
        {
          period=period-1000;
          sum=sum+period*30;
        }
      else if(period>1000)
        {
          period=period-1000;
          sum=sum+period*100;
        }
      else if(period>500)
        {

          period=period-250;
          sum=sum+period*100;
        }
      else if(period>400)
        {
          period=period-50;
          sum=sum+period*100;
        }
      else if(period>300)
        {
          period=period-10;
          sum=sum+period*100;
        }
      else if(period>250)
        {
          period=period-5;
          sum=sum+period*100;
        }
      else if(period>50)
        {
          period=period-1;
          sum=sum+period*100;
        }
      else
        {
          printf("restitution protocol finished\n");
          exit(1);
        }
    }

      if(time>=tbegin && time<=tend)
    {
      Istim=stimstrength;
    }
      if(time>tend)
    {
      Istim=0.;
      tbegin=tbegin+period;
      tend=tbegin+stimduration;
    }
#endif


#ifdef S1S2RESTPROTOCOL
      if(counter<repeats)
    {

      if(time>=tbegin && time<=tend)
        {
          Istim=stimstrength;
        }
      if(time>tend)
        {
          Istim=0.;
          tbegin=tbegin+basicperiod;
          tend=tbegin+stimduration;
          counter++;
        }

    }
      else if(counter==repeats)
    {

      if(time>=tbegin && time<=tend)
        {
          Istim=stimstrength;
        }
      if(time>tend)
        {
          Istim=0.;
          tbegin=tbegin+basicapd+dia;
          tbeginS2=tbegin;
          tend=tbegin+stimduration;
          counter++;
        }
    }
      else if(counter==repeats+1)
    {
      if(time>=tbegin && time<=tend)
        {
          Istim=stimstrength;
        }
      if(time>tend)
        {
          Istim=0.;
          tbegin=tbegin+basicperiod;
          tend=tbegin+stimduration;
          counter=0;


          if(dia>1000)
        {
          dia=dia-1000;
        }
          else if(dia>300)
        {
          dia=dia-100;
        }
          else if(dia>150)
        {
          dia=dia-5;
        }
          else if(dia>5)
        {
          dia=dia-1;
        }
          else
        {
          printf("restitution protocol finished\n");
          exit(1);
        }
        }
    }
#endif


      Step(&Var,HT,despath,&time,step,Istim);

      if(step % 100 ==0)
      {
          Var.writebackup(&time,despath);
      }
    }
  return 0;
}



























void Step(Variables *V,double HT,char *despath,double *tt,int step,double Istim)
{

#define v(array_pointer,i,j) (*(V->array_pointer+i*V->NJ +j))



  double IKr;
  double IKs;
  double IK1;
  double Ito;
  double INa;
  double IbNa;
  double ICaL;
  double IbCa;
  double INaCa;
  double IpCa;
  double IpK;
  double INaK;
  double Irel;
  double Ileak;
  double Iup;
  double Ixfer;
  double k1;
  double k2;
  double kCaSR;


  double dNai;
  double dKi;
  double dCai;
  double dCaSR;
  double dCaSS;
  double dRR;


  double Ek;
  double Ena;
  double Eks;
  double Eca;
  double CaCSQN;
  double bjsr;
  double cjsr;
  double CaSSBuf;
  double bcss;
  double ccss;
  double CaBuf;
  double bc;
  double cc;
  double Ak1;
  double Bk1;
  double rec_iK1;
  double rec_ipK;
  double rec_iNaK;
  double AM;
  double BM;
  double AH_1;
  double BH_1;
  double AH_2;
  double BH_2;
  double AJ_1;
  double BJ_1;
  double AJ_2;
  double BJ_2;
  double M_INF;
  double H_INF;
  double J_INF;
  double TAU_M;
  double TAU_H;
  double TAU_J;
  double axr1;
  double bxr1;
  double axr2;
  double bxr2;
  double Xr1_INF;
  double Xr2_INF;
  double TAU_Xr1;
  double TAU_Xr2;
  double Axs;
  double Bxs;
  double Xs_INF;
  double TAU_Xs;
  double R_INF;
  double TAU_R;
  double S_INF;
  double TAU_S;
  double Ad;
  double Bd;
  double Cd;
  double Af;
  double Bf;
  double Cf;
  double Af2;
  double Bf2;
  double Cf2;
  double TAU_D;
  double D_INF;
  double TAU_F;
  double F_INF;
  double TAU_F2;
  double F2_INF;
  double TAU_FCaSS;
  double FCaSS_INF;


  static double inverseVcF2=1/(2*Vc*F);
  static double inverseVcF=1./(Vc*F);
  static double inversevssF2=1/(2*Vss*F);


  static char s[200];
  FILE *FF;

  // define all variables
#define       sm          (*V).M
#define       sh          (*V).H
#define       sj          (*V).J
#define       sxr1        (*V).Xr1
#define       sxr2        (*V).Xr2
#define       sxs         (*V).Xs
#define       ss          (*V).S
#define       sr          (*V).R
#define       sd          (*V).D
#define       sf          (*V).F
#define       sf2         (*V).F2
#define       sfcass      (*V).FCass
#define       sRR         (*V).RR
#define       sOO         (*V).OO
#define       svolt       (*V).Volt
#define       svolt2      (*V).Volt2
#define       Cai         (*V).Cai
#define       CaSR        (*V).CaSR
#define       CaSS        (*V).CaSS
#define       Nai         (*V).Nai
#define       Ki          (*V).Ki
#define       sItot       (*V).Itot



    //Needed to compute currents
    Ek=RTONF*(log((Ko/Ki)));
    Ena=RTONF*(log((Nao/Nai)));
    Eks=RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));
    Eca=0.5*RTONF*(log((Cao/Cai)));
    Ak1=0.1/(1.+exp(0.06*(svolt-Ek-200)));
    Bk1=(3.*exp(0.0002*(svolt-Ek+100))+
     exp(0.1*(svolt-Ek-10)))/(1.+exp(-0.5*(svolt-Ek)));
    rec_iK1=Ak1/(Ak1+Bk1);
    rec_iNaK=(1./(1.+0.1245*exp(-0.1*svolt*F/(R*T))+0.0353*exp(-svolt*F/(R*T))));
    rec_ipK=1./(1.+exp((25-svolt)/5.98));


    //Compute currents
    INa=GNa*sm*sm*sm*sh*sj*(svolt-Ena);
    ICaL=GCaL*sd*sf*sf2*sfcass*4*(svolt-15)*(F*F/(R*T))*
      (0.25*exp(2*(svolt-15)*F/(R*T))*CaSS-Cao)/(exp(2*(svolt-15)*F/(R*T))-1.);
    Ito=Gto*sr*ss*(svolt-Ek);
    IKr=Gkr*sqrt(Ko/5.4)*sxr1*sxr2*(svolt-Ek);
    IKs=Gks*sxs*sxs*(svolt-Eks);
    IK1=GK1*rec_iK1*(svolt-Ek);
    INaCa=knaca*(1./(KmNai*KmNai*KmNai+Nao*Nao*Nao))*(1./(KmCa+Cao))*
      (1./(1+ksat*exp((n-1)*svolt*F/(R*T))))*
      (exp(n*svolt*F/(R*T))*Nai*Nai*Nai*Cao-
       exp((n-1)*svolt*F/(R*T))*Nao*Nao*Nao*Cai*2.5);
    INaK=knak*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*rec_iNaK;
    IpCa=GpCa*Cai/(KpCa+Cai);
    IpK=GpK*rec_ipK*(svolt-Ek);
    IbNa=GbNa*(svolt-Ena);
    IbCa=GbCa*(svolt-Eca);


    //Determine total current
    (sItot) = IKr    +
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
      IpK   +
      Istim;






    //update concentrations
    kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR)));
    k1=k1_/kCaSR;
    k2=k2_*kCaSR;
    dRR=k4*(1-sRR)-k2*CaSS*sRR;
    sRR+=HT*dRR;
    sOO=k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);


    Irel=Vrel*sOO*(CaSR-CaSS);
    Ileak=Vleak*(CaSR-Cai);
    Iup=Vmaxup/(1.+((Kup*Kup)/(Cai*Cai)));
    Ixfer=Vxfer*(CaSS-Cai);


    CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
    dCaSR=HT*(Iup-Irel-Ileak);
    bjsr=Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
    cjsr=Kbufsr*(CaCSQN+dCaSR+CaSR);
    CaSR=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;


    CaSSBuf=Bufss*CaSS/(CaSS+Kbufss);
    dCaSS=HT*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-ICaL*inversevssF2*CAPACITANCE));
    bcss=Bufss-CaSSBuf-dCaSS-CaSS+Kbufss;
    ccss=Kbufss*(CaSSBuf+dCaSS+CaSS);
    CaSS=(sqrt(bcss*bcss+4*ccss)-bcss)/2;


    CaBuf=Bufc*Cai/(Cai+Kbufc);
    dCai=HT*((-(IbCa+IpCa-2*INaCa)*inverseVcF2*CAPACITANCE)-(Iup-Ileak)*(Vsr/Vc)+Ixfer);
    bc=Bufc-CaBuf-dCai-Cai+Kbufc;
    cc=Kbufc*(CaBuf+dCai+Cai);
    Cai=(sqrt(bc*bc+4*cc)-bc)/2;


    dNai=-(INa+IbNa+3*INaK+3*INaCa)*inverseVcF*CAPACITANCE;
    Nai+=HT*dNai;

    dKi=-(Istim+IK1+Ito+IKr+IKs-2*INaK+IpK)*inverseVcF*CAPACITANCE;
    Ki+=HT*dKi;


    //write currents to file
    if(step%100==0)
      {
    sprintf(s,"%s%s",despath,"/Currents");
    FF=fopen(s,"a");
    fprintf(FF,"%4.10f\t",(*tt));
    fprintf(FF,"%4.10f\t",IKr);
    fprintf(FF,"%4.10f\t",IKs);
    fprintf(FF,"%4.10f\t",IK1);
    fprintf(FF,"%4.10f\t",Ito);
    fprintf(FF,"%4.10f\t",INa);
    fprintf(FF,"%4.10f\t",IbNa);
    fprintf(FF,"%4.10f\t",INaK);
    fprintf(FF,"%4.10f\t",ICaL);
    fprintf(FF,"%4.10f\t",IbCa);
    fprintf(FF,"%4.10f\t",INaCa);
    fprintf(FF,"%4.10f\t",Irel);
    fprintf(FF,"%4.10f\t",Ileak);
    fprintf(FF,"%4.10f",Ixfer);
    fprintf(FF,"\n");
    fclose(FF);
      }

    //compute steady state values and time constants
    AM=1./(1.+exp((-60.-svolt)/5.));
    BM=0.1/(1.+exp((svolt+35.)/5.))+0.10/(1.+exp((svolt-50.)/200.));
    TAU_M=AM*BM;
    M_INF=1./((1.+exp((-56.86-svolt)/9.03))*(1.+exp((-56.86-svolt)/9.03)));
    if (svolt>=-40.)
      {
    AH_1=0.;
    BH_1=(0.77/(0.13*(1.+exp(-(svolt+10.66)/11.1))));
    TAU_H= 1.0/(AH_1+BH_1);
      }
    else
      {
    AH_2=(0.057*exp(-(svolt+80.)/6.8));
    BH_2=(2.7*exp(0.079*svolt)+(3.1e5)*exp(0.3485*svolt));
    TAU_H=1.0/(AH_2+BH_2);
      }
    H_INF=1./((1.+exp((svolt+71.55)/7.43))*(1.+exp((svolt+71.55)/7.43)));
    if(svolt>=-40.)
      {
    AJ_1=0.;
    BJ_1=(0.6*exp((0.057)*svolt)/(1.+exp(-0.1*(svolt+32.))));
    TAU_J= 1.0/(AJ_1+BJ_1);
      }
    else
      {
     AJ_2=(((-2.5428e4)*exp(0.2444*svolt)-(6.948e-6)*
        exp(-0.04391*svolt))*(svolt+37.78)/
           (1.+exp(0.311*(svolt+79.23))));
     BJ_2=(0.02424*exp(-0.01052*svolt)/(1.+exp(-0.1378*(svolt+40.14))));
     TAU_J= 1.0/(AJ_2+BJ_2);
      }
    J_INF=H_INF;

    Xr1_INF=1./(1.+exp((-26.-svolt)/7.));
    axr1=450./(1.+exp((-45.-svolt)/10.));
    bxr1=6./(1.+exp((svolt-(-30.))/11.5));
    TAU_Xr1=axr1*bxr1;
    Xr2_INF=1./(1.+exp((svolt-(-88.))/24.));
    axr2=3./(1.+exp((-60.-svolt)/20.));
    bxr2=1.12/(1.+exp((svolt-60.)/20.));
    TAU_Xr2=axr2*bxr2;

    Xs_INF=1./(1.+exp((-5.-svolt)/14.));
    Axs=(1400./(sqrt(1.+exp((5.-svolt)/6))));
    Bxs=(1./(1.+exp((svolt-35.)/15.)));
    TAU_Xs=Axs*Bxs+80;

#ifdef EPI
    R_INF=1./(1.+exp((20-svolt)/6.));
    S_INF=1./(1.+exp((svolt+20)/5.));
    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
#endif
#ifdef ENDO
    R_INF=1./(1.+exp((20-svolt)/6.));
    S_INF=1./(1.+exp((svolt+28)/5.));
    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    TAU_S=1000.*exp(-(svolt+67)*(svolt+67)/1000.)+8.;
#endif
#ifdef MCELL
    R_INF=1./(1.+exp((20-svolt)/6.));
    S_INF=1./(1.+exp((svolt+20)/5.));
    TAU_R=9.5*exp(-(svolt+40.)*(svolt+40.)/1800.)+0.8;
    TAU_S=85.*exp(-(svolt+45.)*(svolt+45.)/320.)+5./(1.+exp((svolt-20.)/5.))+3.;
#endif


     D_INF=1./(1.+exp((-8-svolt)/7.5));
     Ad=1.4/(1.+exp((-35-svolt)/13))+0.25;
     Bd=1.4/(1.+exp((svolt+5)/5));
     Cd=1./(1.+exp((50-svolt)/20));
     TAU_D=Ad*Bd+Cd;
     F_INF=1./(1.+exp((svolt+20)/7));
     Af=1102.5*exp(-(svolt+27)*(svolt+27)/225);
     Bf=200./(1+exp((13-svolt)/10.));
     Cf=(180./(1+exp((svolt+30)/10)))+20;
     TAU_F=Af+Bf+Cf;
     F2_INF=0.67/(1.+exp((svolt+35)/7))+0.33;
     Af2=600*exp(-(svolt+25)*(svolt+25)/170);
     Bf2=31/(1.+exp((25-svolt)/10));
     Cf2=16/(1.+exp((svolt+30)/10));
     TAU_F2=Af2+Bf2+Cf2;
     FCaSS_INF=0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;
     TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;



     //Update gates
     sm = M_INF-(M_INF-sm)*exp(-HT/TAU_M);
     sh = H_INF-(H_INF-sh)*exp(-HT/TAU_H);
     sj = J_INF-(J_INF-sj)*exp(-HT/TAU_J);
     sxr1 = Xr1_INF-(Xr1_INF-sxr1)*exp(-HT/TAU_Xr1);
     sxr2 = Xr2_INF-(Xr2_INF-sxr2)*exp(-HT/TAU_Xr2);
     sxs = Xs_INF-(Xs_INF-sxs)*exp(-HT/TAU_Xs);
     ss= S_INF-(S_INF-ss)*exp(-HT/TAU_S);
     sr= R_INF-(R_INF-sr)*exp(-HT/TAU_R);
     sd = D_INF-(D_INF-sd)*exp(-HT/TAU_D);
     sf =F_INF-(F_INF-sf)*exp(-HT/TAU_F);
     sf2 =F2_INF-(F2_INF-sf2)*exp(-HT/TAU_F2);
     sfcass =FCaSS_INF-(FCaSS_INF-sfcass)*exp(-HT/TAU_FCaSS);



     //update voltage
     svolt= svolt + HT*(-sItot);


}

