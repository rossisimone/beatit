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
#include "Electrophysiology/IonicModels/Fabbri17.hpp"
#include "Util/IO/io.hpp"
#include <cmath>
#include <iomanip>
#include "Util/CTestUtil.hpp"
#include "Util/IO/io.hpp"
#include "libmesh/getpot.h"


struct Stimulus
{
        double t0;
        double tf;
        bool on;
        double amp;
        double cl;

        Stimulus() : t0(10.0), tf(10.5), on(false), amp(0), cl(1000) {}
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
                if(on) return -amp;
                else return 0.0;
        }
};


int main(int argc, char ** argv)
{
    BeatIt::printBanner(std::cout);
	// Import input file
	GetPot data = BeatIt::readInputFile(argc, argv);

	std::string output_folder = data("output_folder","Output");
	BeatIt::createOutputFolder(output_folder);

	std::unique_ptr<BeatIt::IonicModel> pFabbri17( BeatIt::IonicModel::IonicModelFactory::Create("Fabbri17") );
	int numVar = pFabbri17->numVariables();
	std::vector<double> variables(numVar, 0.0);
	std::vector<double> currents(12, 0.0);

	std::ofstream output(output_folder+"/results.txt");
	std::ofstream output_iion(output_folder+"/currents.txt");
	pFabbri17->initializeSaveData(output);

	pFabbri17->setup(data,"san");
	pFabbri17->initialize(variables);

//	double C = data("C", 5.7e-5);
//	pFabbri17->set_membrane_capacitance(C);

        std::cout << std::setprecision(15);
        // model is in seconds
	double dt = 0.001; // = 1ms
        dt =1e-2;
	int save_iter = static_cast<int>(1.0 / dt);
//	save_iter = 1;
        double TF = 10e3; //500;//2*1000;
        //TF = 1e-2;
	double Ist = 0;
	double time = 0.0;
	BeatIt::saveData(0.0, variables, output);

	int iter = 0;
	Stimulus stimulus;
	// for ctest purposes
	double solution_norm = 0.0;
	while( time <= TF )
	{
		Ist = stimulus.get(time);

		pFabbri17->solve(variables, Ist, dt);
		time += dt;
		++iter;
		if( 0 == iter%save_iter )
		{
			BeatIt::saveData(time, variables, output);
			pFabbri17->get_currents(currents);;
			BeatIt::saveData(time, currents, output_iion);
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

