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
 * \brief This test provides the port of the Grandi 2011 model from Matlab to c++
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
#include <memory>
#include "Electrophysiology/IonicModels/Grandi11.hpp"
#include "Util/IO/io.hpp"
#include "Util/CTestUtil.hpp"
#include <iomanip>

int main()
{
    BeatIt::printBanner(std::cout);

    std::unique_ptr<BeatIt::IonicModel> pORd( BeatIt::IonicModel::IonicModelFactory::Create("Grandi11") );
    int numVar = pORd->numVariables();
    std::vector<double> variables(numVar, 0.0);

    std::ofstream output("results.txt");
    pORd->initializeSaveData(output);
    pORd->initialize(variables);

    double dt = 6e-3;

    int save_iter = 1; //1 / dt;
    double TF = 300.0;
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
            Ist = 80.0;
        }
        else Ist = 0.0;

        pORd->solve(variables, Ist, dt);
        time += dt;
        ++iter;
        if( 0 == iter%save_iter ) BeatIt::saveData(time, variables, output);

        // for ctest purposes
        solution_norm += variables[0];

    }
    output.close();
    //for ctest purposes
    solution_norm /= iter;
    std::cout << std::setprecision(18) << "Solution norm = " << solution_norm << std::endl;
    //up to the 16th digit
    const double reference_solution_norm = -48.0678940914975072;
    //We check only up to 12th
	return BeatIt::CTest::check_test(solution_norm, reference_solution_norm, 1e-12);

}
