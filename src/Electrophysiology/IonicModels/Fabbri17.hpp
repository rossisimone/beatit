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
 * \file Fabbri17.hpp
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

#ifndef SRC_ELECTROPHYSIOLOGY_IONICMODELS_FABBRI17_HPP_
#define SRC_ELECTROPHYSIOLOGY_IONICMODELS_FABBRI17_HPP_

#include "Electrophysiology/IonicModels/IonicModel.hpp"

namespace BeatIt
{

/*!
 *
 */
class Fabbri17: public IonicModel
{
public:
	/// Ionic Model
	typedef IonicModel super;
	//! Constructor of the Fabbri 2017 SAN model
	/*!
	 *
	 */
	Fabbri17();


	//! Update all the variables in the ionic model
	/*!
	 *  \param [in] variables Vector containing the local value of all variables
	 *  \param [in] dt        Timestep
	 */
	void updateVariables(std::vector<double>& variables, double appliedCurrent,
			double dt);
	//! Update all the variables in the ionic model
	/*!
	 *  \param [in] V transmember potential (Variables  does not include potential)
	 *  \param [in] variables Vector containing the local value of all variables
	 *  \param [in] dt        Timestep
	 */
	void updateVariables(double V, std::vector<double>& variables, double dt) {}

	//! Evaluate total ionic current for the computation of the potential
	/*!
	 *  \param [in] variables Vector containing the local value of all variables
	 *  \param [in] appliedCurrent value of the applied current
	 *  \param [in] dt        Timestep
	 */
	double evaluateIonicCurrent(std::vector<double>& variables,
			double appliedCurrent = 0.0, double dt = 0.0);
	//! Evaluate total ionic current for the computation of the potential
	/*!
	 *  \param [in] V transmember potential
	 *  \param [in] variables Vector containing the local value of all variables (Variables  does not  include potential)
	 *  \param [in] appliedCurrent value of the applied current
	 *  \param [in] dt        Timestep
	 */
	double evaluateIonicCurrent(double V, std::vector<double>& variables,
			double appliedCurrent = 0.0, double dt = 0.0)
	{
		return 0.0;
	}

    // NOT YET IMPLEMENTED
	double evaluateIonicCurrentTimeDerivative(std::vector<double>& variables,
			std::vector<double>& old_variables, double dt = 0.0,
			double h = 0.0);

	//! Initialize the values of the variables
	/*!
	 *  \param [in] variables Vector containing the local value of all variables
	 */
	void initialize(std::vector<double>& variables);

	//! Initialize the output files with the names
	/*!
	 *  \param [in] output output file where we will save the values
	 */
	void initializeSaveData(std::ostream& output);
	//Nernst Potentials
	void nernst_potentials(const std::vector<double>& variables);

private:
	void update();
	/* Ion Current Functions */


	/* Cell Geometry */
    double C; /* cell capacitance (pF) */
    double L_cell;// 67 /* Length of the cell (um) */
    constexpr static double L_sub = 0.02; /* Distance between jSR and membrane surface (submembrane space) (um) */
    constexpr static double pi = 3.14159265358979323846; /* Pi */
    double R_cell; // 3.9 /* Radius of the cell (um) */
    constexpr static double V_ipart =  0.46; /* part of the cell volume occupied by myoplasm */
    constexpr static double V_jsrpart =  0.0012; /* part of the cell volume occupied by junctional SR */
    constexpr static double V_nsrpart =  0.46; /* part of the cell volume occupied by network SR */
    double V_cell;// = pi * R_cell * R_cell * L_cell; /*cell volume */
    double V_sub;// = 2 * pi * L_sub * L_cell * ( R_cell - 0.5 * L_sub ); /* submembrane space volume */
    double V_i;// = V_ipart * V_cell - V_sub; /*myoplasmic  volume */
    double V_jsr;// = V_jsrpart * V_cell; /*JSR  volume */
    double V_nsr;// = V_nsrpart * V_cell; /*NSR  volume */
    /* done */

    /* Fixed ion concentrations, mM. */
    constexpr static double Cao = 1.8; // extracellular Ca2+ concentration
    constexpr static double Ki = 140.0;  // intracellular K+ concentration
    constexpr static double Ko = 5.4;  // extracellular K+ concentration
    constexpr static double Nao = 140.0; // extracellular Na+ concentration
    //constexpr static double Mgi = 5.0; // intracellular Mg2+ concentration
    /* done */

    double E_Na; //reversal potential for Na+
    double E_mh; //reversal potential for fast Na+ channel
    double E_K; //reversal potential for K+
    double E_Ks;//reversal potential for slow rectifier K+ channel


    /* Terms for Solution of Conductance and Reversal Potential */
    constexpr static double Rc = 8314.472; /* Universal Gas Constant (J/kmol*K) */
    constexpr static double F = 96485.0; /* Faraday's Constant (C/mol) */
    constexpr static double T = 310.0; /* Temperature (K) */
    constexpr static double RTONF = Rc * T / F; // (mV)
    /* done */


    // Sarcolemmal ion currents and their conductances
    // If
    const double g_f = 0.00427; //uS
//    const double g_fNa = 0.00268; //uS
//    const double g_fK = 0.00159; //uS
    //ICaL
    const double P_CaL = 0.4578; // nA/mM
    //ICaT
    const double P_CaT = 0.04132; // nA/mM
    // IKr
    const double g_Kr = 0.00424; //uS
    // IKs
    const double g_Ks_ = 0.00065; //uS  // gks 0.000299 XXXXXXXXXXXX
    // IKACh
    const double g_KACh = 0.00345; // uS
    //Ito
    const double g_to = 0.0035; //uS
    // INa
    const double g_Na = 0.0223; // uS
    const double g_Na_L = 0; // uS
    // INaCa
    const double K_NaCa = 3.343; // nA
    // IKur
    const double g_Kur = 0.0001539; // uS
    // SONO ARRIVATO QUI Jan 03 2019!

    // Modulation of sarcolemmal ion currents by ions
    const double Km_fCa = 0.000338; // mM - dissiociation constant of Ca2+ dependent ICal inactivation
    const double Km_Kp = 1.4; // mM - half maximal Ko for INaK v
    const double Km_Nap = 14.0; // mM - half maximal Nai for INaK v
    const double alpha_fCa = 0.0075; // #define alpha_fca 0.021 // s^-1 Ca2+ dissociation rate constant for ICal
    const double i_NaK_max = 0.08105; // nanoA
    // Na+/Ca2+ exchangere (NaCa) function
    const double K1ni = 395.3; //mM v
    const double K1no = 1628.0; //mM v
    const double K2ni = 2.289; //mM v
    const double K2no = 561.4; //mM v
    const double K3ni = 26.44; //mM v
    const double K3no = 4.663; //mM v
    const double Kci = 0.0207; //mM v
    const double Kcni = 26.44; //mM v
    const double Kco = 3.663; //mM v
    const double Qci = 0.1369; //mM v
    const double Qco = 0.0; //mM v
    const double Qn = 0.4315; //mM v

    // Ca2+ diffusion
    const double taudifCa =  5.469e-5; // s
    const double tautr = 0.04; //s

    // SERCA pump
//    const double Kup = 0.000286113;  286.; // mM // Kup 0.0006 XXXXXXXXXXXXXXXXX
  //  const double Pup = 5.; // mM/s
    //const double slopeup = 50.0; //mM

    // RyR function
    const double kiCa = 500.0; // 1/mM/s
    const double kim = 5.0; // 1/s  // kim 0.005
    const double koCa = 10000.0; // 1/mM^2 /s
    const double kom = 660; // 1/s
    const double ks = 1.48e8; // 1/s
    const double EC50_SR = 0.45; //mM // v
    const double HSR = 2.5; // v
    const double MaxSR = 15.; // v
    const double MinSR = 1.;  // v

    // Ca2+ and Mg2+ buffering
    const double CMtot = 0.045; // mM
    const double CQtot = 10; // mM
    const double TCtot = 0.031; // mM
    const double TMCtot = 0.062; // mM
    const double kbCM = 542; // 1/s
    const double kbCQ = 445; // 1/s
    const double kbTC = 446; // 1/s
    const double kbTMC = 7.51; // 1/s
    const double kbTMM = 751; // 1/s
    const double kfCM = 1.642e6; // 1/s /mM
    const double kfCQ = 175.4; // 1/s /mM
    const double kfTC = 88800; // 1/s /mM
    const double kfTMC = 227700; // 1/s/mM
    const double kfTMM = 2277; // 1/s/mM

    // Funny current
    const double Km_f =  45; //mM
    const double alpha =  0.5927;
    const double blockade = 0;


	double v; /* Membrane voltage (mV) */
	double vnew; /* New Voltage (mV) */
	double dvdt; /* Change in Voltage / Change in Time (mV/ms) */
	double dvdtnew; /* New dv/dt (mV/ms) */
	double dt;      /* Time step (ms) */
	double t;       /* Time (ms) */
	double istim;       /* Constant Stimulus (uA/cm^2) */
	double itot;       /* Total current (uA/cm^2) */

    double i_f;// nanoA {pub; //in};
    double i_NaK; //nanoA {pub; //in};
    double i_NaCa; //nanoA {pub; //in};
    double i_Na; //nanoA {pub; //in};
    double i_NaL; //nanoA {pub; //in};
    double i_Kr; //nanoA {pub; //in};
    double i_Ks; //nanoA {pub; //in};
    double i_to; //nanoA {pub; //in};
    double i_CaL; //nanoA {pub; //in};
    double i_CaT; //nanoA {pub; //in};
    double i_KACh; //nanoA {pub; //in};
    double i_Kur; //nanoA {pub; //in};
    double V; //millivolt {pub; //out};
    double clamp_mode; //dimensionless {init; //0};
    double V_clamp; //millivolt {priv; //in};
    double V_ode; //millivolt {init; //-47.787168};
    double i_tot; //nanoA;

    double Nai_clamp;
    double i_fNa; // nanoA {pub: in};
    double i_siNa; // nanoA {pub: in};
    double dnai;

    double blockade_NaCa;
    double doo;
    double di;

    double j_SRCarel;
    double P_tot;

    const double slope_up = 5e-5;
    const double P_up_basal = 5;
    const double tau_dif_Ca = 5.469e-5;
    const double tau_tr = 0.04;
    const double K_up = 0.000286113;

    // Ca_buffering
    double TC_tot = 0.031; // ConcTC
    double TMC_tot = 0.062; // ConcTMC
    double CM_tot = 0.045; // ConcCM 0.045
    double CQ_tot = 10.0; // ConcCQ 10.0
    double kf_TC = 88800; // kfTC 88.8
    double kf_TMM = 2277; // kfTMM 2.277
    double kf_TMC = 227700; // kfTMC 237.7
    double kf_CM = 1.642e6; // kfCM 237.7  XXXXXXXXXXXXXXXX
    double kf_CQ = 175.4;  // kfCQ 0.534  XXXXXXXXXXXXXXXXXX
    double kb_TC = 446.0; // kbTC 0.446
    double kb_TMC = 7.51; // kbTMC 0.00751
    double kb_TMM = 751.0; // kbTMM 0.751
    double kb_CM = 542.0; // kbCM 0.542
    double kb_CQ = 445.0;  // kbCQ 0.445
    const double Mgi = 2.5; // v

    double k_dL =  4.3371;
    double V_dL = -16.4508;

    // Nai
    // 1
    double Nai;// = 5.0; // intracellular Na+ concentration


    //funny
    // 2
    double y; // 0.009508
    double y_infinity;
    double y_shift;
    //Ina
    // 3
    double m; //0.447724
    const double delta_m = 1e-5;
    // 4
    double h; //0.003058
    //CaL
    // 5
    double dL; //0.001921
    // 6
    double fL; //0.846702
    // 7
    double fCa; // 0.844449

    // CaT
    // 8
    double dT;//0.268909
    // 9
    double fT;//0.020484

    //Ca_SR_release
    // 10
    double R; // 0.9308
    // 11
    double O; // 6.181512e-9
    // 12
    double I; // 4.595622e-10
    // 13
    double RI; // 0.069199

    // Ca buffering
    // 14
    double fTC; //0.017929
    // 15
    double fTMC; //0.259947
    // 16
    double fTMM; //0.653777
    // 17
    double fCMi; //0.217311
    // 18
    double fCMs; //0.158521
    // 19
    double fCQ; //0.138975

    // Ca dynamics
    // 20
    double Cai; //9.15641e-6
    // 21
    double Ca_sub; // 6.226104e-5
    // 22
    double Ca_nsr; //0.435148
    // 23
    double Ca_jsr; // 0.409551

    // Ikur
    // 24
    double r_Kur; // 0.011845
    // 25
    double s_Kur; // 0.845304

    // i_to
    // 26
    double q; // 0.430836
    // 27
    double r; // 0.014523

    // I_Kr
    // 28
    double paS; // 0.283185
    // 29
    double paF; // 0.011068
    // 30
    double piy; // 0.709051

    // I_KS
    // 31
    double n; // 0.1162

    // I_KACh a gate
    // 32
    double a; //0.00277

    void get_currents(std::vector<double>& currents);

};


/*! Creator function: it calls the constructor of ORd
 *
 *  @return IonicModel* pointer to a new ORd object
 */
IonicModel* createFabbri17();

namespace
{
static bool register_Fabbri17 = IonicModel::IonicModelFactory::Register(
		"Fabbri17", &createFabbri17);
}

} /* namespace BeatIt */

#endif /* SRC_ELECTROPHYSIOLOGY_IONICMODELS_Fabbri17_HPP_ */
