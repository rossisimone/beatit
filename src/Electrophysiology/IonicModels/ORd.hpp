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
 * \file ORd.hpp
 *
 * \class ORd
 *
 * \brief This class provides the implementationa of the OHara-Rudy model
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

#ifndef SRC_IONICMODELS_ORD_HPP_
#define SRC_IONICMODELS_ORD_HPP_

#include "Electrophysiology/IonicModels/IonicModel.hpp"

namespace BeatIt
{

class ORd: public IonicModel
{
public:
	/// Ionic Model
	typedef IonicModel 						super;
	//! Constructor of the OHara-Rudy Ionic Model
	/*!
	 *
	 */
	ORd();
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
    void updateVariables(std::vector<double>& variables, std::vector<double>& rhs, double appliedCurrent, double dt, bool overwrite);
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
private:


    /// original methods from OHara Rudy code
    void revpots(std::vector<double>& variables);
    void RGC(std::vector<double>& variables, double dt);
    void FBC(std::vector<double>& variables, double dt);

    /// constants
    constexpr const static double nao = 140.0;//extracellular sodium in mM
    constexpr const static double cao = 1.8;//extracellular calcium in mM
    constexpr const static double ko  = 5.4;//extracellular potassium in mM

    /// buffer paramaters
    constexpr const static double BSRmax  = 0.047;
    constexpr const static double KmBSR   = 0.00087;
    constexpr const static double BSLmax  = 1.124;
    constexpr const static double KmBSL   = 0.0087;
    constexpr const static double cmdnmax = 0.05;
    constexpr const static double kmcmdn  = 0.00238;
    constexpr const static double trpnmax = 0.07;
    constexpr const static double kmtrpn  = 0.0005;
    constexpr const static double csqnmax = 10.0;
    constexpr const static double kmcsqn  = 0.8;

    /// CaMK paramaters
    constexpr const static double aCaMK  = 0.05;
    constexpr const static double bCaMK  = 0.00068;
    constexpr const static double CaMKo  = 0.05;
    constexpr const static double KmCaM  = 0.0015;
    constexpr const static double KmCaMK = 0.15;

    /// cell geometry
    constexpr const static double L     = 0.01;
    constexpr const static double rad   = 0.0011;
    constexpr const static double vcell = 1000*3.14*rad*rad*L;
    constexpr const static double Ageo  = 2*3.14*rad*rad+2*3.14*rad*L;
    constexpr const static double Acap  = 2*Ageo;
    constexpr const static double vmyo  = 0.68*vcell;
    constexpr const static double vmito = 0.26*vcell;
    constexpr const static double vsr   = 0.06*vcell;
    constexpr const static double vnsr  = 0.0552*vcell;
    constexpr const static double vjsr  = 0.0048*vcell;
    constexpr const static double vss   = 0.02*vcell;


    /// introduce varaibles for reversal potentials, currents, fluxes, and CaMK
    double ENa,EK,EKs;
    double INa,INaL,Ito,ICaL,ICaNa,ICaK,IKr,IKs,IK1,INaCa_i,INaCa_ss,INaCa,INaK,IKb,INab,IpCa,ICab,Ist;
    double Jrel,Jup,Jtr,Jdiff,JdiffNa,JdiffK,Jleak;
    double CaMKa,CaMKb;
};


/*! Creator function: it calls the constructor of ORd
 *
 *  @return IonicModel* pointer to a new ORd object
 */
IonicModel* createORd();

namespace
{
    static bool register_ORd = IonicModel::IonicModelFactory::Register("ORd", &createORd );
}


} /* namespace BeatIt */

#endif /* SRC_IONICMODELS_ORD_HPP_ */
