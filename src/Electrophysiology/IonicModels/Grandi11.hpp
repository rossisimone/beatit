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
 * \file Grandi11.hpp
 *
 * \class Grandi11
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
class Grandi11: public IonicModel
{
public:
    /// Ionic Model
    typedef IonicModel                      super;
    //! Constructor of the OHara-Rudy Ionic Model
    /*!
     *
     */
    Grandi11();
    //! Solve method
    /*!
     *  \param [in] variables Vector containing the local value of all variables
     *  \param [in] appliedCurrent value of the applied current
     *  \param [in] dt        Timestep
     */
    void solve(std::vector<double>& variables, double appliedCurrent, double dt);
    void setup(GetPot& data, std::string section);

    //! Update all the variables in the ionic model
    /*!
     *  \param [in] variables Vector containing the local value of all variables
     *  \param [in] dt        Timestep
     */
    void updateVariables(std::vector<double>& variables, double appliedCurrent, double dt);
    //! Update all the variables in the ionic model
    /*!
     *  \param [in] V transmember potential (Variables  does not include potential)
     *  \param [in] variables Vector containing the local value of all variables
     *  \param [in] dt        Timestep
     */
   void updateVariables(double V, std::vector<double>& variables, double dt){}


    //! Evaluate total ionic current for the computation of the potential
    /*!
     *  \param [in] variables Vector containing the local value of all variables
     *  \param [in] appliedCurrent value of the applied current
     *  \param [in] dt        Timestep
     */
    double evaluateIonicCurrent(std::vector<double>& variables, double appliedCurrent = 0.0, double dt = 0.0);
    //! Evaluate total ionic current for the computation of the potential
    /*!
     *  \param [in] V transmember potential
     *  \param [in] variables Vector containing the local value of all variables (Variables  does not  include potential)
     *  \param [in] appliedCurrent value of the applied current
     *  \param [in] dt        Timestep
     */
     double evaluateIonicCurrent(double V, std::vector<double>& variables, double appliedCurrent = 0.0, double dt = 0.0){ return 0.0;}
     double evaluatedIonicCurrent(std::vector<double>& variables, double appliedCurrent = 0.0, double dt = 0.0, double h = 0.0);
     double evaluatedIonicCurrent(std::vector<double>& variables, std::vector<double>& old_variables, double dt = 0.0, double h = 0.0);

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
    /// Physical Constants
    constexpr static double S_PI = 3.14159265358979323846264338327950288;
    constexpr static double Frdy = F;
    constexpr static double Temp = T;
    constexpr static double FoRT = Frdy/R/Temp;
    // should it be 1.381?
    constexpr static double Cmem = 1.1e-10; // 1.3810e-10;  // membrane capacitance 1.3810e-10 [F/m^2];%
    constexpr static double Qpow = (Temp-310.0)/10.0;

    /// Cell geometry
    constexpr static double cellLength = 100;     // cell length [um]113;//100
    constexpr static double cellRadius = 10.25;   // cell radius [um]12;//10.25
    constexpr static double junctionLength = 160e-3;  // junc length [um]
    constexpr static double junctionRadius = 15e-3;   // junc radius [um]
    constexpr static double distSLcyto = 0.45;    // dist. SL to cytosol [um]
    constexpr static double distJuncSL = 0.5;  // dist. junc to SL [um]
    constexpr static double DcaJuncSL = 1.64e-6;  // Dca junc to SL [cm^2/sec]
    constexpr static double DcaSLcyto = 1.22e-6; // Dca SL to cyto [cm^2/sec]
    constexpr static double DnaJuncSL = 1.09e-5;  // Dna junc to SL [cm^2/sec]
    constexpr static double DnaSLcyto = 1.79e-5;  // Dna SL to cyto [cm^2/sec]
    constexpr static double Vcell = S_PI*cellRadius*cellRadius*cellLength*1e-15;    // [L]
    constexpr static double Vmyo = 0.65*Vcell;
    constexpr static double Vsr = 0.035*Vcell;
    constexpr static double Vsl = 0.02*Vcell;
    constexpr static double Vjunc = 1*0.0539*.01*Vcell;
    constexpr static double SAjunc = 20150*S_PI*2*junctionLength*junctionRadius;  // [um^2]
    constexpr static double SAsl = S_PI*2*cellRadius*cellLength;          // [um^2]

    ///
    constexpr static double J_ca_juncsl = 1. / 1.2134e12; // [L/msec] = 8.2413e-13
    constexpr static double J_ca_slmyo  = 1. / 2.68510e11; // [L/msec] = 3.2743e-12
    constexpr static double J_na_juncsl = 1. / (1.6382e12/3.*100); // [L/msec] = 6.1043e-13
    constexpr static double J_na_slmyo  = 1. / (1.8308e10/3.*100);  // [L/msec] = 5.4621e-11

    //  Fractional currents in compartments
    constexpr static double Fjunc = 0.11;
    constexpr static double Fsl = 1-Fjunc;
    constexpr static double Fjunc_CaL = 0.9;
    constexpr static double Fsl_CaL = 1-Fjunc_CaL;

    // Fixed ion concentrations
    constexpr static double Cli = 15;   // Intracellular Cl  [mM]
    constexpr static double Clo = 150;  // Extracellular Cl  [mM]
    constexpr static double Ko = 5.4;   // Extracellular K   [mM]
    constexpr static double Nao = 140;  // Extracellular Na  [mM]
    constexpr static double Cao = 1.8;  // Extracellular Ca  [mM]
    constexpr static double Mgi = 1;    // Intracellular Mg  [mM]

    //Nernst Potentials
    double ena_junc, ena_sl, ek, eca_junc, eca_sl, ecl;

    //
    constexpr static double AF = 0.0;
    constexpr static double ISO = 0.0;

    // Na transport parameters
    constexpr static double GNa=23*(1-0.1*AF);  // [mS/uF]
    constexpr static double GNaB = 0.597e-3;   // [mS/uF]
    constexpr static double IbarNaK = 1.26;     // [uA/uF]
    constexpr static double KmNaip = 11*(1-0.25*ISO);         // [mM]11
    constexpr static double KmKo =1.5;         // [mM]1.5
    constexpr static double Q10NaK = 1.63;
    constexpr static double Q10KmNai = 1.39;

    // K current parameters
    constexpr static double pNaK = 0.01833;
    constexpr static double gkp = 0.002;

    // Cl current parameters
    constexpr static double GClCa =0.0548;   // [mS/uF]
    constexpr static double GClB = 9e-3;        // [mS/uF]
    constexpr static double KdClCa = 100e-3;    // [mM]

    // I_Ca parameters
    constexpr static double pNa = (1+0.5*ISO)*(1-0.5*AF)*0.75e-8;       // [cm/sec]
    constexpr static double pCa = (1+0.5*ISO)*(1-0.5*AF)*2.7e-4;       // [cm/sec]
    constexpr static double pK = (1+0.5*ISO)*(1-0.5*AF)*1.35e-7;        // [cm/sec]
    constexpr static double Q10CaL = 1.8;


    // Ca transport parameters
    constexpr static double IbarNCX = (1+0.4*AF)*3.15;      // [uA/uF]5.5 before - 9 in rabbit
    constexpr static double KmCai = 3.59e-3;    // [mM]
    constexpr static double KmCao = 1.3;        // [mM]
    constexpr static double KmNai = 12.29;      // [mM]
    constexpr static double KmNao = 87.5;       // [mM]
    constexpr static double ksat = 0.27;        // [none]
    constexpr static double nu = 0.35;          // [none]
    constexpr static double Kdact =0.384e-3;   // [mM] 0.256 rabbit
    constexpr static double Q10NCX = 1.57;      // [none]
    constexpr static double IbarSLCaP =  0.0471; // IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
    constexpr static double KmPCa =0.5e-3;     // [mM]
    constexpr static double GCaB = 6.0643e-4;    // [uA/uF] 3
    constexpr static double Q10SLCaP = 2.35;    // [none]

    // SR flux parameters
    constexpr static double Q10SRCaP = 2.6;          // [none]
    constexpr static double Vmax_SRCaP = 5.3114e-3;  // [mM/msec] (286 umol/L cytosol/sec)
    constexpr static double Kmf = (2.5-1.25*ISO)*0.246e-3;          // [mM] default
    constexpr static double Kmr = 1.7;               // [mM]L cytosol
    constexpr static double hillSRCaP = 1.787;       // [mM]
    constexpr static double ks = 25;                 // [1/ms]
    constexpr static double koCa = 10+20*AF+10*ISO*(1-AF);               // [mM^-2 1/ms]   //default 10   modified 20
    constexpr static double kom = 0.06;              // [1/ms]
    constexpr static double kiCa = 0.5;              // [1/mM/ms]
    constexpr static double kim = 0.005;             // [1/ms]
    constexpr static double ec50SR = 0.45;           // [mM]

    // Buffering parameters
    // koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
    constexpr static double Bmax_Naj = 7.561;       // [mM] // Na buffering
    constexpr static double Bmax_Nasl = 1.65;       // [mM]
    constexpr static double koff_na = 1e-3;         // [1/ms]
    constexpr static double kon_na = 0.1e-3;        // [1/mM/ms]
    constexpr static double Bmax_TnClow = 70e-3;    // [mM]                      // TnC low affinity
    constexpr static double koff_tncl = (1+0.5*ISO)*19.6e-3;    // [1/ms]
    constexpr static double kon_tncl = 32.7;        // [1/mM/ms]
    constexpr static double Bmax_TnChigh = 140e-3;  // [mM]                      // TnC high affinity
    constexpr static double koff_tnchca = 0.032e-3; // [1/ms]
    constexpr static double kon_tnchca = 2.37;      // [1/mM/ms]
    constexpr static double koff_tnchmg = 3.33e-3;  // [1/ms]
    constexpr static double kon_tnchmg = 3e-3;      // [1/mM/ms]
    constexpr static double Bmax_CaM = 24e-3;       // [mM] **? about setting to 0 in c-code**   // CaM buffering
    constexpr static double koff_cam = 238e-3;      // [1/ms]
    constexpr static double kon_cam = 34;           // [1/mM/ms]
    constexpr static double Bmax_myosin = 140e-3;   // [mM]                      // Myosin buffering
    constexpr static double koff_myoca = 0.46e-3;   // [1/ms]
    constexpr static double kon_myoca = 13.8;       // [1/mM/ms]
    constexpr static double koff_myomg = 0.057e-3;  // [1/ms]
    constexpr static double kon_myomg = 0.0157;     // [1/mM/ms]
    constexpr static double Bmax_SR = 19*.9e-3;     // [mM] (Bers text says 47e-3) 19e-3
    constexpr static double koff_sr = 60e-3;        // [1/ms]
    constexpr static double kon_sr = 100;           // [1/mM/ms]
    constexpr static double Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        // [mM]    // SL buffering
    constexpr static double Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    // [mM]    //Fei *0.1!!! junction reduction factor
    constexpr static double koff_sll = 1300e-3;     // [1/ms]
    constexpr static double kon_sll = 100;          // [1/mM/ms]
    constexpr static double Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       // [mM]
    constexpr static double Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  // [mM] //Fei *0.1!!! junction reduction factor
    constexpr static double koff_slh = 30e-3;       // [1/ms]
    constexpr static double kon_slh = 100;          // [1/mM/ms]
    constexpr static double Bmax_Csqn = 140e-3*Vmyo/Vsr;            // [mM] // Bmax_Csqn = 2.6;      // Csqn buffering
    constexpr static double koff_csqn = 65;         // [1/ms]
    constexpr static double kon_csqn = 100;         // [1/mM/ms]
    constexpr static double GNaL=0.0025*AF;


    //
    double I_Na, I_Na_junc, I_Na_sl;
    double I_NaL, I_NaL_junc, I_NaL_sl;
    double I_nabk, I_nabk_junc, I_nabk_sl;
    double I_nak, I_nak_junc, I_nak_sl;
    double I_kr;
    double I_ks, I_ks_junc, I_ks_sl;
    double I_kp, I_kp_junc, I_kp_sl;

    //I_to
    // %nS/pF maleckar; %human atrium
    constexpr static double GtoFast = 0.165*(1.0-0.7*AF);
    double I_to, I_tof;

    double I_kur, I_ki;
    double I_ClCa_junc, I_ClCa_sl, I_ClCa;
    double I_Clbk, I_ClCFTR;
    double I_Ca_junc, I_Ca_sl, I_Ca;
    double I_CaK;
    double I_CaNa_junc, I_CaNa_sl, I_CaNa;
    double I_Catot;
    double I_ncx_junc, I_ncx_sl, I_ncx;
    double I_pca_junc, I_pca_sl, I_pca;
    double I_cabk_junc, I_cabk_sl, I_cabk;

    double J_SRCarel, J_serca, J_SRleak;
    double J_CaB_cytosol, J_CaB_junction, J_CaB_sl;

    double I_Na_tot_junc, I_Na_tot_sl, I_Na_tot_sl2, I_Na_tot_junc2;
    double I_K_tot;
    double I_Ca_tot_junc, I_Ca_tot_sl;

    double I_Na_tot, I_Cl_tot, I_Ca_tot;
    double I_tot;
    bool set_resting_conditions;
    bool markov_iks;
};


/*! Creator function: it calls the constructor of ORd
 *
 *  @return IonicModel* pointer to a new ORd object
 */
IonicModel* createGrandi11();

namespace
{
    static bool register_Grandi11 = IonicModel::IonicModelFactory::Register("Grandi11", &createGrandi11 );
}



} /* namespace BeatIt */

#endif /* SRC_ELECTROPHYSIOLOGY_IONICMODELS_GRANDI11_HPP_ */
