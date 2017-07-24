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
 * \file main_ORd.cpp
 *
 * \class main_ORd
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
 * Created on: Aug 7, 2016
 *
 */

#include <memory>
#include "Electrophysiology/IonicModels/TP06.hpp"
#include "Util/IO/io.hpp"
#include <cmath>
#include <iomanip>
#include "Util/CTestUtil.hpp"
#include "Electromechanics/Activation/NegroniLascano.hpp"
#include "libmesh/getpot.h"

struct Stimulus
{
        double t0;
        double tf;
        bool on;
        double amp;
        double cl;

        Stimulus() : t0(0.0), tf(0.5), on(false), amp(80), cl(500) {}
        Stimulus(double T0, double TF,  double cycleLength) : t0(T0), tf(T0+TF), on(false), amp(80), cl(cycleLength) {}
        double get(double time)
        {
                if (time < t0)
                        on = false;
                else if (time >= t0 && time <= tf)
                        on = true;
                else // time > tf
                {
                        t0 += cl;
                        tf += cl;
                        on = false;
                }
                if(on) return amp;
                else return 0.0;
        }
};


int main( int argc, char ** argv)
{
    BeatIt::printBanner(std::cout);
    auto data = BeatIt::readInputFile(argc, argv);

    double t0 = data("t0", 0.0);
    double duration = data("duration", 0.5);
    double tf = data("tf", 1000.0);
    double cl = data("cl", 1000.0);
	double dt =data("dt", 0.00025);
	std::string ionic_model = data("model", "NONE");
	std::string output_name = data("output", "results");
	double save_dt = data("saveDt", 0.00025);

	std::unique_ptr<BeatIt::IonicModel> pORd( BeatIt::IonicModel::IonicModelFactory::Create(ionic_model) );
	int numVar = pORd->numVariables();
	std::vector<double> variables(numVar, 0.0);
	std::ofstream output(output_name+".txt");
	pORd->initializeSaveData(output);
	pORd->initialize(variables);


	BeatIt::NegroniLascano nl;
	std::ofstream output_nl(output_name+"_nl.txt");
	int nlVar = nl.numVariables();
	std::vector<double> nl_variables(nlVar, 0.0);
	nl.initializeSaveData(output_nl);

	int save_iter = save_dt / dt;
    double TF = tf;
	double Ist = 0;
	double time = 0.0;
	BeatIt::saveData(0.0, variables, output);
	BeatIt::saveData(0.0, nl_variables, output_nl);

	int iter = 0;
	Stimulus stimulus(t0, duration, cl);
	// for ctest purposes
	double solution_norm = 0.0;
	while( time <= TF )
	{
		Ist = stimulus.get(time);
		pORd->solve(variables, Ist, dt);
		double Cai = 0.0;
		if("TP06" ==  ionic_model) Cai = 1e3*variables[1];
		else if("Grandi11" ==  ionic_model) Cai = 1e3*variables[36];
		double I4f = 1.1;
		if(time > 30) I4f = 1.1 - 0.25 * time / 1000;
		double Cai_diast = 0.0;
		nl.solve(nl_variables, Cai, dt, I4f, Cai_diast);
		time += dt;
		++iter;
		if( 0 == iter%save_iter )
		{
			BeatIt::saveData(time, variables, output);
			BeatIt::saveData(time, nl_variables, output_nl);
		}

		// for ctest purposes
		solution_norm += variables[0];

	}
	output.close();
	//for ctest purposes
	solution_norm /= iter;
	//up to the 16th digit
	std::cout << std::setprecision(18) << "Solution norm = " << solution_norm << std::endl;
	const double reference_solution_norm = -18.2035079050909516;
	//We check only up to 12th
	return BeatIt::CTest::check_test(solution_norm, reference_solution_norm, 1e-10);
}

