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
 * \file FentonKarma.hpp
 *
 * \class FentonKarma
 *
 * \brief This class provides the implementation of the Nash Panfilov model
 *
 * For details on thi model refer to:
 * Nash, M.P. and Panfilov, A.V., 2004.
 * Electromechanical model of excitable tissue to study reentrant cardiac arrhythmias.
 * Progress in biophysics and molecular biology, 85(2), pp.501-522.
 *
 *
 * \author srossi
 *
 * \version 0.0
 *
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Aug 11, 2016
 *
 */

#ifndef SRC_ELECTROPHYSIOLOGY_IONICMODELS_FentonKarma_HPP_
#define SRC_ELECTROPHYSIOLOGY_IONICMODELS_FentonKarma_HPP_

#include "Electrophysiology/IonicModels/IonicModel.hpp"

namespace BeatIt
{

/*!
 *
 */
class FentonKarma: public IonicModel
{
public:
    /// Ionic Model
    typedef IonicModel                      super;
    //! Constructor of the OHara-Rudy Ionic Model
    /*!
     *
     */
    FentonKarma();
    ~FentonKarma() {};

    //! Update all the variables in the ionic model
    /*!
     *  \param [in] variables Vector containing the local value of all variables (Variables  includes potential)
     *  \param [in] dt        Timestep
     */
    void updateVariables(std::vector<double>& variables, double appliedCurrent, double dt);
    void updateVariables(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent, double dt) ;
    void updateVariables(std::vector<double>& variables, std::vector<double>& rhs, double appliedCurrent, double dt, bool overwrite);
    bool isSecondOrderImplemented() { return true; }
    //! Update all the variables in the ionic model
    /*!
     *  \param [in] V transmember potential (Variables  does not include potential)
     *  \param [in] variables Vector containing the local value of all variables
     *  \param [in] dt        Timestep
     */
   void updateVariables(double V, std::vector<double>& variables, double dt);


    //! Evaluate total ionic current for the computation of the potential
    /*!
     *  \param [in] variables Vector containing the local value of all variables (Variables  includes potential)
     *  \param [in] appliedCurrent value of the applied current
     *  \param [in] dt        Timestep
     */
    double evaluateIonicCurrent(std::vector<double>& variables, double appliedCurrent, double dt);
    double evaluateIonicCurrent(std::vector<double>& v_n, std::vector<double>& v_np1, double appliedCurrent = 0.0, double dt = 0.0);
	//! Evaluate total ionic current for the computation of the potential
	/*!
     *  \param [in] V transmember potential
	 *  \param [in] variables Vector containing the local value of all variables (Variables  does not  include potential)
	 *  \param [in] appliedCurrent value of the applied current
	 *  \param [in] dt        Timestep
	 */
     double evaluateIonicCurrent(double V, std::vector<double>& variables, double appliedCurrent = 0.0, double dt = 0.0);

     double evaluateIonicCurrentTimeDerivative(std::vector<double>& variables, std::vector<double>& old_variables, double dt = 0.0, double h = 0.0);


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

    void setup(GetPot& data, std::string section = "monodomain");

    double evaluateSAC(double v , double I4f);
private:

    double M_tau_v_p;
    double M_tau_v1_m;
    double M_tau_v2_m;
    double M_tau_w_p;
    double M_tau_w_m;
    double M_tau_d;
    double M_tau_0;
    double M_tau_r;
    double M_tau_si;
    double M_kappa;
    double M_V_c_si;
    double M_V_c;
    double M_V_v;

    int p(double V) { return (V >= M_V_c ) ?  1 : 0; }
    int q(double V) { return (V >= M_V_v ) ?  1 : 0; }
};


/*! Creator function: it calls the constructor of ORd
 *
 *  @return IonicModel* pointer to a new ORd object
 */
IonicModel* createFentonKarma();

namespace
{
    static bool register_FentonKarma = IonicModel::IonicModelFactory::Register("FentonKarma", &createFentonKarma );
}



} /* namespace BeatIt */

#endif /* SRC_ELECTROPHYSIOLOGY_IONICMODELS_FentonKarma_HPP_ */
