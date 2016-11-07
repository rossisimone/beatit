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


struct Pacing
{
    double istim(double t)
    {
        update(t);
        if(on) return amplitude;
        else return 0.0;
    }

    void update(double t)
    {
        if(t > end_time && on == true)
        {
            on = false;
            start_time += cycle_length;
            end_time += cycle_length;
        }
        if(t >= start_time && t <= end_time)
        {
            on = true;
        }


    }
    bool on;
    double start_time;
    double end_time;
    double cycle_length;
    double amplitude;

};


void checkAPD(bool& on, double APD_threshold, double time, std::vector<double>& APD, double& APD_time, double v)
{
		if(on)
		{
			if( v < APD_threshold)
			{
				APD_time -= time;
//				std::cout << "Finishing an AP at time: " << time << ", with APD: " << APD_time << std::endl;
				APD.push_back(-APD_time);
				APD_time = 0.0;
				on = false;
			}
		}
		else
		{
			if( v >= APD_threshold)
			{
//				std::cout << "Starting an AP at time: " << time << std::endl;
				on = true;
				APD_time = time;
			}
		}
}

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
	double TF = 3000.0;
	double Ist = 0;
	double time = 0.0;
	BeatIt::saveData(0.0, variables, output);

	int iter = 0;


	// for ctest purposes
	double solution_norm = 0.0;

	Pacing pacing;
	pacing.amplitude = 1.0;
    pacing.start_time= 0.0;
    pacing.end_time= 0.5;
    pacing.cycle_length= 15.0;

    bool on75 = false;
    bool on50 = false;
    bool on25 = false;
    std::vector<double> APD75;
    std::vector<double> APD50;
    std::vector<double> APD25;
    double APD75_time = 0.0;
    double APD50_time = 0.0;
    double APD25_time = 0.0;
    double APD75_threshold = 0.75;
    double APD50_threshold = 0.5;
    double APD25_threshold = 0.25;

    while( time <= TF )
	{
	    Ist = pacing.istim(time);
//		if(time >= 0. && time <= 0.5)
//		{
//			Ist = 1.0;
//		}
//		else Ist = 0.0;

		pNP->solve(variables, Ist, dt);
		time += dt;
		++iter;
		if( 0 == iter%save_iter ) BeatIt::saveData(time, variables, output);

		// for ctest purposes
		solution_norm += variables[0];

		checkAPD(on75, APD75_threshold, time, APD75, APD75_time, variables[0]);
		checkAPD(on50, APD50_threshold, time, APD50, APD50_time, variables[0]);
		checkAPD(on25, APD25_threshold, time, APD25, APD25_time, variables[0]);
	}

	output.close();
	//for ctest purposes
	solution_norm /= iter;
	//up to the 16th digit
	const double reference_solution_norm =    0.6704271748046145;

	std::cout << "APD75 size: " << APD75.size() << std::endl;
//	for(auto&& v : APD75) std::cout << v << std::endl;
		std::cout << APD75[APD75.size()-1] << std::endl;
	std::cout << "APD50 size: " << APD50.size() << std::endl;
//	for(auto&& v : APD50) std::cout << v << std::endl;
		std::cout << APD50[APD50.size()-1] << std::endl;
	std::cout << "APD25 size: " << APD25.size() << std::endl;
//	for(auto&& v : APD25)
		std::cout << APD25[APD25.size()-1] << std::endl;
//	for(int i = 0; i < APD.size(); ++i)
//	{
//		std::cout << "APD in cycle: " << i+1 << APD[i] << std::endl;
//	}
	return BeatIt::CTest::check_test(solution_norm, reference_solution_norm, 1e-12);

}

