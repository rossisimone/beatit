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
 * \file TP06.hpp
 *
 * \class TP06
 *
 * \brief This class provides the implementationa of the 10Tusscher-Panfilov 2006 model
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

#ifndef SRC_IONICMODELS_TP06_HPP_
#define SRC_IONICMODELS_TP06_HPP_

#include "Electrophysiology/IonicModels/IonicModel.hpp"

namespace BeatIt
{

class TP06: public IonicModel
{
public:
	/// Ionic Model
	typedef IonicModel 						super;
	//! Constructor of the OHara-Rudy Ionic Model
	/*!
	 *
	 */
	TP06();
    //! Solve method
	/*!
	 *  \param [in] variables Vector containing the local value of all variables
	 *  \param [in] appliedCurrent value of the applied current
	 *  \param [in] dt        Timestep
	 */
    void solve(std::vector<double>& variables, double appliedCurrent, double dt);

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

    void step(std::vector<double>& variables, double dt);
    void setCellType(CellType type);
    void selectParameters(CellType type);

    double membraneCapacitance()
    {
        return 1.0;
    }
private:



    /*-----------------------------------------------------------------------------
      ELECTROPHYSIOLOGICAL PARAMETERS:
    -----------------------------------------------------------------------------*/

    //External concentrations
    constexpr static double Ko =5.4;
    constexpr static double Cao=2.0;
    constexpr static double Nao=140.0;

    //Intracellular volumes
    constexpr static double Vc =0.016404;
    constexpr static double Vsr=0.001094;
    constexpr static double Vss=0.00005468;

    constexpr static double inverseVcF2 =1. / (2 * Vc * F);
    constexpr static double inverseVcF  =1. / (Vc * F);
    constexpr static double inversevssF2=1. / (2. * Vss * F);

    //Calcium buffering dynamics
    constexpr static double Bufc  =0.2;
    constexpr static double Kbufc =0.001;
    constexpr static double Bufsr =10.;
    constexpr static double Kbufsr=0.3;
    constexpr static double Bufss =0.4;
    constexpr static double Kbufss=0.00025;

    //Intracellular calcium flux dynamics
    constexpr static double Vmaxup=0.006375;
    constexpr static double Kup   =0.00025;
    constexpr static double Vrel  =0.102;//40.8;
    constexpr static double k1_   =0.15;
    constexpr static double k2_   =0.045;
    constexpr static double k3    =0.060;
    constexpr static double k4    =0.005;//0.000015;
    constexpr static double EC    =1.5;
    constexpr static double maxsr =2.5;
    constexpr static double minsr =1.;
    constexpr static double Vleak =0.00036;
    constexpr static double Vxfer =0.0038;

    //Cellular capacitance
    constexpr static double CAPACITANCE=0.185;

    //Parameters for currents
    //Parameters for IKr
    constexpr static double Gkr=0.153;
    //Parameters for Iks
    constexpr static double pKNa=0.03;

    double Gks;
    double Gto;
    //Parameters for Ik1
    constexpr static double GK1=5.405;

    //Parameters for INa
    constexpr static double GNa=14.838;
    //Parameters for IbNa
    constexpr static double GbNa=0.00029;
    //Parameters for INaK
    constexpr static double KmK=1.0;
    constexpr static double KmNa=40.0;
    constexpr static double knak=2.724;
    //Parameters for ICaL
    constexpr static double GCaL=0.00003980;
    //Parameters for IbCa
    constexpr static double GbCa=0.000592;
    //Parameters for INaCa
    constexpr static double knaca=1000;
    constexpr static double KmNai=87.5;
    constexpr static double KmCa=1.38;
    constexpr static double ksat=0.1;
    constexpr static double n=0.35;
    //Parameters for IpCa
    constexpr static double GpCa=0.1238;
    constexpr static double KpCa=0.0005;
    //Parameters for IpK;
    constexpr static double GpK=0.0146;

    double Itot;
    double dItot;
    double Istim;
    double Volt2;

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
};


/*! Creator function: it calls the constructor of TP06
 *
 *  @return IonicModel* pointer to a new TP06 object
 */
IonicModel* createTP06();

namespace
{
    static bool register_TP06 = IonicModel::IonicModelFactory::Register("TP06", &createTP06 );
}


} /* namespace BeatIt */

#endif /* SRC_IONICMODELS_TP06_HPP_ */
