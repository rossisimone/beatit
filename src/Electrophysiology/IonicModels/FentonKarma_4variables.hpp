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
 * Comparison of Detailed and Simplified Models of Human Atrial Myocytes to Recapitulate Patient Specific Properties
 * Daniel M. Lombardo, Flavio H. Fenton, Sanjiv M. Narayan, Wouter-Jan Rappel
 * Computational Biology
 *
 * NOTE that this model is somewhat different than the original
 *  Pulmonary vein reentry—Properties and size matter: Insights from a computational analysis
 *
 * \author srossi
 *
 * \version 0.0
 *
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Mar 02, 2021
 *
 */

#ifndef SRC_ELECTROPHYSIOLOGY_IONICMODELS_FentonKarma_4variables_HPP_
#define SRC_ELECTROPHYSIOLOGY_IONICMODELS_FentonKarma_4variables_HPP_

#include "Electrophysiology/IonicModels/IonicModel.hpp"

namespace BeatIt
{

/*!
 *
 */
class FentonKarma4v: public IonicModel
{
public:
    /// Ionic Model
    typedef IonicModel                      super;
    //! Constructor of the OHara-Rudy Ionic Model
    /*!
     *
     */
    FentonKarma4v();
    ~FentonKarma4v() {};

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

    double M_tv1m; // ms
    double M_tv2m; // ms
    double M_tvp; // ms

    double M_tw1m; // ms
    double M_tw2m; // ms
    double M_twp; // ms

    double M_to;  // ms
    double M_ta;  // ms
    double M_td;  // ms
    double M_tsi; // ms


    double M_uc;
    double M_uv;
    double M_uw;
    double M_uo;
    double M_um;

    double M_ucsi;
    double M_uso;

    double M_rsm;
    double M_rsp;

    double M_k;

    double M_aso;// ms
    double M_bso;// ms
    double M_cso;// ms

    bool M_pulmonary_vein_model;

    double H(double x) { return (x > 0 ) ?  1.0 : 0.0; }
};


/*! Creator function: it calls the constructor of ORd
 *
 *  @return IonicModel* pointer to a new ORd object
 */
IonicModel* createFentonKarma4v();

namespace
{
    static bool register_FentonKarma4v = IonicModel::IonicModelFactory::Register("FentonKarma4v", &createFentonKarma4v );
}



} /* namespace BeatIt */

#endif /* SRC_ELECTROPHYSIOLOGY_IONICMODELS_FentonKarma_HPP_ */
