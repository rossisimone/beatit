/*
 * Activation.hpp
 *
 *  Created on: Jul 13, 2017
 *      Author: srossi
 */

#ifndef SRC_ELECTROMECHANICS_ACTIVATION_ACTIVATION_HPP_
#define SRC_ELECTROMECHANICS_ACTIVATION_ACTIVATION_HPP_

#include <vector>
#include <string>
#include "Util/Factory.hpp"


namespace BeatIt
{

class Activation
{
public:
public:

	/// Create a factory
	typedef Factory<Activation, std::string>     ActivationFactory;


	Activation(int numVar, const std::string name = "empty");
	virtual ~Activation();

    virtual void solve(std::vector<double>& variables, double Cai, double dt, double I4f = 1.0, double Cai_diast = 0.0) = 0;
	//! Initialize the values of the variables
	/*!
	 *  \param [in] variables Vector containing the local value of all variables
	 */
	virtual void initialize(std::vector<double>& variables) = 0;

	//! Initialize the output files with the names
	/*!
	 *  \param [in] output output file where we will save the values
	 */
    virtual void initializeSaveData(std::ostream& output) = 0;

    //virtual void setup(GetPot& data, std::string section) {}
    int numVariables()
     {
     	return M_numVariables;
     }


    /// Total number of variables, including the potential
    int    M_numVariables;
    /// Name of the variables (excluding the potential)
    std::vector<std::string> M_variablesNames;
    /// Name of the ionic model
    std::string M_modelName;
};

} /* namespace BeatIt */

#endif /* SRC_ELECTROMECHANICS_ACTIVATION_ACTIVATION_HPP_ */
