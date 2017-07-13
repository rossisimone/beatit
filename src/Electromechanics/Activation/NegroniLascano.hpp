/*
 * NegroniLascano.hpp
 *
 *  Created on: Jul 13, 2017
 *      Author: srossi
 */

#ifndef SRC_ELECTROMECHANICS_ACTIVATION_NEGRONILASCANO_HPP_
#define SRC_ELECTROMECHANICS_ACTIVATION_NEGRONILASCANO_HPP_

#include "Electromechanics/Activation/Activation.hpp"

namespace BeatIt
{

class NegroniLascano: public Activation
{
public:
	typedef Activation super;
	NegroniLascano();
	virtual ~NegroniLascano();

    void solve(std::vector<double>& variables, double Cai, double dt, double I4f = 1.0, double Cai_diast = 0.0);
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

    //void setup(GetPot& data, std::string section) {}


	constexpr static double Y1 = 39 ; // uM/s
	constexpr static double Z1 = 30 ; // /s
	constexpr static double Y2 = 1.3 ; // /s
	constexpr static double Z2 = 1.3 ; // /s
	constexpr static double Y3 = 30 ; // /s
	constexpr static double Z3 = 1560 ; // uM/s
	constexpr static double Y4 = 40 ; // /s
	constexpr static double Yd = 9 ; // s/uM^2
	constexpr static double Tt = 70 ; // uM
	constexpr static double B = 1200 ; // /s
	constexpr static double hc = 0.005 ; // um
	constexpr static double La = 1.17 ; // um
	constexpr static double R = 20 ; // /um^2
	constexpr static double L0 = 0.97 ; // um


	constexpr static double A = 944.5815 ; // mN/mm^2/um/uM (adjusted)
    //A = 0.0625*A ; // mN/mm^2/um/uM (adjusted)
	constexpr static double B_s = 0.4853 ; // mN/mm^2 (adjusted)
	constexpr static double h = 3.0843/B_s ; // cm (Kass89)

    double TCa_eff;
    double T, TCa, TCas, Ts, X, F;
    double dTCa, dTCas, dTs, dX;
    double Qb, Qa, Qr, Qd, Qd1, Qd2;
};

} /* namespace BeatIt */

#endif /* SRC_ELECTROMECHANICS_ACTIVATION_NEGRONILASCANO_HPP_ */
