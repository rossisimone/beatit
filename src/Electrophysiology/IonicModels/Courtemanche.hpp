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
 * \file Courtemanche.hpp
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

#ifndef SRC_ELECTROPHYSIOLOGY_IONICMODELS_GRANDI11_HPP_
#define SRC_ELECTROPHYSIOLOGY_IONICMODELS_GRANDI11_HPP_

#include "Electrophysiology/IonicModels/IonicModel.hpp"

namespace BeatIt
{

/*!
 *
 */
class Courtemanche: public IonicModel
{
public:
	/// Ionic Model
	typedef IonicModel super;
	//! Constructor of the OHara-Rudy Ionic Model
	/*!
	 *
	 */
	Courtemanche();
	//! Solve method
	/*!
	 *  \param [in] variables Vector containing the local value of all variables
	 *  \param [in] appliedCurrent value of the applied current
	 *  \param [in] dt        Timestep
	 */
	void solve(std::vector<double>& variables, double appliedCurrent,
			double dt);

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
	double evaluatedIonicCurrent(std::vector<double>& variables,
			double appliedCurrent = 0.0, double dt = 0.0, double h = 0.0) { return 0.0; }
	double evaluatedIonicCurrent(std::vector<double>& variables,
			std::vector<double>& old_variables, double dt = 0.0,
			double h = 0.0) { return 0.0; }

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
	void comp_ina(); /* Calculates Fast Na Current */
	void comp_ical(); /* Calculates Currents through L-Type Ca Channel */
	void comp_ikr(); /* Calculates Rapidly Activating K Current */
	void comp_iks(); /* Calculates Slowly Activating K Current */
	void comp_iki(); /* Calculates Time-Independant K Current */
	void comp_ikach(); /* Calculates Acetylcholine-sensitive potassium*/
	void comp_ikur(); /* Calculates Ultra-Rapidly activation K Current*/
	void comp_ito(); /* Calculates Transient Outward Current */
	void comp_inaca(); /* Calculates Na-Ca Exchanger Current */
	void comp_inak(); /* Calculates Na-K Pump Current */
	void comp_ipca(); /* Calculates Sarcolemmal Ca Pump Current */
	void comp_icab(); /* Calculates Ca Background Current */
	void comp_inab(); /* Calculates Na Background Current */
	void comp_it(); /* Calculates Total Current */

	/* Ion Concentration Functions */
	void conc_nai(); /* Calculates new myoplasmic Na ion concentration */
	void conc_ki(); /* Calculates new myoplasmic K ion concentration */
	void conc_nsr(); /* Calculates new NSR Ca ion concentration */
	void conc_jsr(); /* Calculates new JSR Ca ion concentration */
	void calc_itr(); /* Calculates Translocation of Ca from NSR to JSR */
	void conc_cai(); /* Calculates new myoplasmic Ca ion concentration */

	/* Cell Geometry */
	constexpr static double l = 0.01; /* Length of the cell (cm) */
	constexpr static double a = 0.0008; /* Radius of the cell (cm) */
	constexpr static double pi = 3.141592; /* Pi */
	/* Cell volume (uL) */
	constexpr static double vcell = 1000 * pi * a * a * l; /*   3.801e-5 uL */
	/* Geometric membrane area (cm^2) */
	constexpr static double ageo = 2 * pi * a * a + 2 * pi * a * l;
	/* Capacity */
	constexpr static double acap = ageo * 2; /*   1.534e-4 cm^2 */
	/* Myoplasm volume (uL) */
	constexpr static double vmyo = vcell * 0.68;
	/* Mitochondria volume (uL) */
	constexpr static double vmito = vcell * 0.26;
	/* SR volume (uL) */
	constexpr static double vsr = vcell * 0.06;
	/* NSR volume (uL) */
	constexpr static double vnsr = vcell * 0.0552;
	/* JSR volume (uL) */
	constexpr static double vjsr = vcell * 0.0048;

	double v; /* Membrane voltage (mV) */
	double vnew; /* New Voltage (mV) */
	double dvdt; /* Change in Voltage / Change in Time (mV/ms) */
	double dvdtnew; /* New dv/dt (mV/ms) */
	double dt;      /* Time step (ms) */
	double t;       /* Time (ms) */
	double st;       /* Constant Stimulus (uA/cm^2) */
	double it;       /* Total current (uA/cm^2) */

	/* Terms for Solution of Conductance and Reversal Potential */
	const double R = 8314; /* Universal Gas Constant (J/kmol*K) */
	const double frdy = 96485; /* Faraday's Constant (C/mol) */
	const double temp = 310; /* Temperature (K) */

	/* Ion Valences */
	const double zna = 1; /* Na valence */
	const double zk = 1; /* K valence */
	const double zca = 2; /* Ca valence */

	/* Ion Concentrations */
	double nai; /* Intracellular Na Concentration (mM) */
	constexpr static double nao = 140; /* Extracellular Na Concentration (mM) */
	double ki; /* Intracellular K Concentration (mM) */
	constexpr static double ko = 4.5; /* Extracellular K Concentration (mM) */
	double cai; /* Intracellular Ca Concentration (mM) */
	constexpr static double cao = 1.8; /* Extracellular Ca Concentration (mM) */
	double cmdn; /* Calmodulin Buffered Ca Concentration (mM) */
	double trpn; /* Troponin Buffered Ca Concentration (mM) */
	double nsr; /* NSR Ca Concentration (mM) */
	double jsr; /* JSR Ca Concentration (mM) */
	double csqn; /* Calsequestrin Buffered Ca Concentration (mM) */

	/* Myoplasmic Na Ion Concentration Changes */
	double naiont; /* Total Na Ion Flow (mM/ms) */
	double dnai; /* Change in Intracellular Na Concentration (mM) */

	/* Myoplasmic K Ion Concentration Changes */
	double kiont; /* Total K Ion Flow (mM/ms) */
	double dki; /* Change in Intracellular K Concentration (mM) */

	/* NSR Ca Ion Concentration Changes */
	double dnsr; /* Change in [Ca] in the NSR (mM) */
	double iup; /* Ca uptake from myo. to NSR (mM/ms) */
	double ileak; /* Ca leakage from NSR to myo. (mM/ms) */
	double kleak; /* Rate constant of Ca leakage from NSR (ms^-1) */
	constexpr static double kmup = 0.00092; /* Half-saturation concentration of iup (mM) */
	constexpr static double iupbar = 0.005; /* Max. current through iup channel (mM/ms) */
	constexpr static double nsrbar = 15; /* Max. [Ca] in NSR (mM) */

	/* JSR Ca Ion Concentration Changes */
	double djsr; /* Change in [Ca] in the JSR (mM) */
	double urel; /* Activation gate u of Ca release from jsr*/
	double urelss; /* Steady state of activation gate u*/
	double tauurel; /* Time constant of activation gate u*/
	double vrel; /* Activation gate v of Ca release from jsr*/
	double vrelss; /* Steady state of activation gate v*/
	double tauvrel; /* Time constant of activation gate v*/
	double wrel; /* Inactivation gate w of Ca release from jsr*/
	double wrelss; /* Steady state of inactivation gate w*/
	double tauwrel; /* Time constant of inactivation gate w*/
	double fn;
	const double grelbarjsrol = 30; /* Rate constant of Ca release from JSR due to overload (ms^-1)*/
	double greljsrol; /* Rate constant of Ca release from JSR due to CICR (ms^-1)*/
	double ireljsrol; /* Ca release from JSR to myo. due to JSR overload (mM/ms)*/
	const double csqnbar = 10; /* Max. [Ca] buffered in CSQN (mM)*/
	const double kmcsqn = 0.8; /* Equalibrium constant of buffering for CSQN (mM)*/
	double bjsr; /* b Variable for analytical computation of [Ca] in JSR (mM)*/
	double cjsr; /* c Variable for analytical computation of [Ca] in JSR (mM)*/
	double on; /* Time constant of activation of Ca release from JSR (ms)*/
	double off; /* Time constant of deactivation of Ca release from JSR (ms)*/
	double magrel; /* Magnitude of Ca release*/

	/* Translocation of Ca Ions from NSR to JSR */
	double itr; /* Translocation current of Ca ions from NSR to JSR (mM/ms)*/
	const double tautr = 180; /* Time constant of Ca transfer from NSR to JSR (ms)*/

	/* Myoplasmic Ca Ion Concentration Changes */
	double caiont; /* Total Ca Ion Flow (mM/ms) */
	double dcai; /* Change in myoplasmic Ca concentration (mM) */
	double b1cai;
	double b2cai;
	constexpr static double cmdnbar = 0.050; /* Max. [Ca] buffered in CMDN (mM) */
	constexpr static double trpnbar = 0.070; /* Max. [Ca] buffered in TRPN (mM) */
	constexpr static double kmcmdn = 0.00238; /* Equalibrium constant of buffering for CMDN (mM) */
	constexpr static double kmtrpn = 0.0005; /* Equalibrium constant of buffering for TRPN (mM) */

	/* Fast Sodium Current (time dependant) */
	double ina; /* Fast Na Current (uA/uF) */
	double gna; /* Max. Conductance of the Na Channel (mS/uF) */
	double ena; /* Reversal Potential of Na (mV) */
	double ah; /* Na alpha-h rate constant (ms^-1) */
	double bh; /* Na beta-h rate constant (ms^-1) */
	double aj; /* Na alpha-j rate constant (ms^-1) */
	double bj; /* Na beta-j rate constant (ms^-1) */
	double am; /* Na alpha-m rate constant (ms^-1) */
	double bm; /* Na beta-m rate constant (ms^-1) */
	double h; /* Na activation */
	double j; /* Na inactivation */
	double m; /* Na inactivation */
	double gB;

	/* Current through L-type Ca Channel */
	double ilca; /* Ca current through L-type Ca channel (uA/uF) */
	double ilcatot; /* Total current through the L-type Ca channel (uA/uF) */
	double ibarca; /* Max. Ca current through Ca channel (uA/uF) */
	double d; /* Voltage dependant activation gate */
	double dss; /* Steady-state value of activation gate d  */
	double taud; /* Time constant of gate d (ms^-1) */
	double f; /* Voltage dependant inactivation gate */
	double fss; /* Steady-state value of inactivation gate f */
	double tauf; /* Time constant of gate f (ms^-1) */
	double fca; /* Ca dependant inactivation gate */
	double taufca; /* Time constant of gate fca (ms^-1) */
	double fcass; /* Steady-state value of activation gate fca  */

	constexpr static double gcalbar = 0.1238;

	/* Acetylcholine-Activated Potassium Current */
	/* modified from Matsuoka et al., Jap J Physiol 2003;53:105-123 */
	double ikach; /* Acetylcholine-activated K current (uA/uF) */
	double gkach; /* Channel conductance of acetylcholine-activated K current (mS/uF) */
	double ekach; /* Reversal potential of acetylcholine-activated K current (mV) */
	double alphayach; /* Alpha rate constant (ms^-1) */
	double betayach; /* Beta rate constant (ms^-1) */
	double tauyach; /* Time constant (ms) */
	double yachss; /* Steady-state value */
	double yach;
	constexpr static double ach = 0.0; /* Acetylcholine concentration */

	/* Ultra-Rapidly Activating Potassium Current */
	double ikur; /* Ultra-rapidly activating K current (uA/uF) */
	double gkur; /* Channel conductance of ultra-rapidly activating K current (mS/uF) */
	double ekur; /* Reversal potential of ultra-rapidly activating K current (mV) */
	double uakur; /* Ultra-rapidly activating K activation gate ua */
	double uakurss; /* Steady-state value of activation gate ua */
	double tauuakur; /* Time constant of gate ua (ms^-1) */
	double alphauakur; /* Alpha rate constant of activation gate ua (ms^-1) */
	double betauakur; /* Beta rate constant of activation gate ua (ms^-1) */
	double uikur; /* Ultra-rapidly activating K activation gate ui*/
	double uikurss; /* Steady-state value of activation gate ui */
	double tauuikur; /* Time constant of gate ui (ms) */
	double alphauikur; /* Alpha rate constant of activation gate ui (ms^-1) */
	double betauikur; /* Beta rate constant of activation gate ui (ms^-1) */

	/* Rapidly Activating Potassium Current */
	double ikr; /* Rapidly activating K current (uA/uF) */
	double gkr; /* Channel conductance of rapidly activating K current (mS/uF) */
	double ekr; /* Reversal potential of rapidly activating K current (mV) */
	double xr; /* Rapidly activating K time-dependant activation */
	double xrss; /* Steady-state value of inactivation gate xr */
	double tauxr; /* Time constant of gate xr (ms^-1) */
	double r; /* K time-independant inactivation */

	/* Slowly Activating Potassium Current */
	double iks; /* Slowly activating K current (uA/uF) */
	double gks; /* Channel conductance of slowly activating K current (mS/uF) */
	double eks; /* Reversal potential of slowly activating K current (mV) */
	double xs; /* Slowly activating potassium current activation gate*/
	double xsss; /* Steady-state value of activation gate xs */
	double tauxs; /* Time constant of gate xs (ms^-1) */
	constexpr static double prnak = 0.01833; /* Na/K Permiability Ratio */

	/* Time-Independent Potassium Current */
	/*Partly modified from Matsuoka, et al, Jap J Physiol,2003:53:105-123*/
	double iki; /* Time-independant K current (uA/uF) */
	double gki; /* Channel conductance of time independant K current (mS/uF) */
	double eki; /* Reversal potential of time independant K current (mV) */
	double kin; /* K inactivation */
	double iku; /*Attaching rate constant of Magnesium to iki*/
	double ikl; /*Detaching rate constant of Magnesium to iki*/
	double ikay; /*Attaching rate constant of spermine to iki*/
	double ikby; /*Detaching rate constant of spermine to iki*/
	double tauiky; /*Time constant of spermine attachment*/
	double ikyss; /*Steady state of spermine attachment*/
	double iky; /*Spermine attachment*/
	double foiki; /*Fraction of channel free from attachment of Magnesium*/
	double fbiki; /*Fraction of channel with attachment of Magnesium*/

	/* Transient Outward Potassium Current */
	double ito; /* Transient outward current */
	double gito; /* Maximum conductance of Ito */
	double erevto; /* Reversal potential of Ito */
	double ato; /* Ito activation */
	double alphaato; /* Ito alpha-a rate constant */
	double betaato; /* Ito beta-a rate constant */
	double tauato; /* Time constant of a gate */
	double atoss; /* Steady-state value of a gate */
	double iito; /* Ito inactivation */
	double alphaiito; /* Ito alpha-i rate constant */
	double betaiito; /* Ito beta-i rate constant */
	double tauiito; /* Time constant of i gate */
	double iitoss; /* Steady-state value of i gate */

	/* Sodium-Calcium Exchanger */
	double inaca; /* NaCa exchanger current (uA/uF) */
	constexpr static double kmnancx = 87.5; /* Na saturation constant for NaCa exchanger */
	constexpr static double ksatncx = 0.1; /* Saturation factor for NaCa exchanger */
	constexpr static double kmcancx = 1.38; /* Ca saturation factor for NaCa exchanger */
	constexpr static double gammas = 0.35; /* Position of energy barrier controlling voltage dependance of inaca */

	/* Sodium-Potassium Pump */
	double inak; /* NaK pump current (uA/uF) */
	double fnak; /* Voltage-dependance parameter of inak */
	double sigma; /* [Na]o dependance factor of fnak */
	constexpr static double ibarnak = 1.0933; /* Max. current through Na-K pump (uA/uF) */
	constexpr static double kmnai = 10; /* Half-saturation concentration of NaK pump (mM) */
	constexpr static double kmko = 1.5; /* Half-saturation concentration of NaK pump (mM) */

	/* Sarcolemmal Ca Pump */
	double ipca; /* Sarcolemmal Ca pump current (uA/uF) */
	constexpr static double ibarpca = 0.275; /* Max. Ca current through sarcolemmal Ca pump (uA/uF) */
	constexpr static double kmpca = 0.0005; /* Half-saturation concentration of sarcolemmal Ca pump (mM) */

	/* Ca Background Current */
	double icab; /* Ca background current (uA/uF) */
	double gcab; /* Max. conductance of Ca background (mS/uF) */
	double ecan; /* Nernst potential for Ca (mV) */

	/* Na Background Current */
	double inab; /* Na background current (uA/uF) */
	double gnab; /* Max. conductance of Na background (mS/uF) */
	double enan; /* Nernst potential for Na (mV) */

	/* Total Ca current */
	double icatot;
	//    double I_tot;
};

/*! Creator function: it calls the constructor of ORd
 *
 *  @return IonicModel* pointer to a new ORd object
 */
IonicModel* createCourtemanche();

namespace
{
static bool register_Courtemanche = IonicModel::IonicModelFactory::Register(
		"Courtemanche", &createCourtemanche);
}

} /* namespace BeatIt */

#endif /* SRC_ELECTROPHYSIOLOGY_IONICMODELS_GRANDI11_HPP_ */
