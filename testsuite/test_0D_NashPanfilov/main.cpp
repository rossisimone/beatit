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
 * \file main.cpp
 *
 * \class main
 *
 * \brief This test check the Nash Panfilov Ionic Model
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

#include <memory>
#include "Electrophysiology/IonicModels/NashPanfilov.hpp"
#include "Util/IO/io.hpp"
#include "Util/CTestUtil.hpp"

int main()
{
    BeatIt::printBanner(std::cout);

	std::unique_ptr<BeatIt::IonicModel> pNP( BeatIt::IonicModel::IonicModelFactory::Create("NashPanfilov") );
	int numVar = pNP->numVariables();
	std::vector<double> variables(numVar, 0.0);

	std::ofstream output("results.txt");
	pNP->initializeSaveData(output);
	pNP->initialize(variables);

	double dt = 0.005;
	int save_iter = 0.1 / dt;
	double TF = 50.0;
	double Ist = 0;
	double time = 0.0;
	BeatIt::saveData(0.0, variables, output);

	int iter = 0;

	// for ctest purposes
	double solution_norm = 0.0;

	while( time <= TF )
	{
		if(time >= 0. && time <= 0.5)
		{
			Ist = 1.0;
		}
		else Ist = 0.0;

		pNP->solve(variables, Ist, dt);
		time += dt;
		++iter;
		if( 0 == iter%save_iter ) BeatIt::saveData(time, variables, output);

		// for ctest purposes
		solution_norm += variables[0];

	}

	output.close();
	//for ctest purposes
	solution_norm /= iter;
	//up to the 16th digit
	const double reference_solution_norm =    0.6704271748046145;

	return BeatIt::CTest::check_test(solution_norm, reference_solution_norm, 1e-12);

}

