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
#include <iomanip>
#include "libmesh/getpot.h"
/*
 * We start by defining an object that will help us keeping track of the pacing protocol
 * In this case we use a simple S1 protocol, that is we apply the same stimulus every time
 * and the time between two stimuli is always the same.
 */
struct Pacing
{
    /*! Get the magnitude of the stimulus at the current time
     *
     * @param t current time
     * @return amplitue of the stimulus
     */
    double istim(double t)
    {
        // First update the Pacing object
        // See description of the update method for more info
        update(t);
        // If stimulus is active the return it's amplitude
        if(on) return amplitude;
        // else return 0
        else return 0.0;
    }

    /*!  Get the current time and update the members of this class
     *
     * @param t current time
     */
    void update(double t)
    {
        // If the stimulus is on
        // but we passed the stimulus final time
        // we need to turn it off
        // In addition we update start and end times
        // for the next stimulus
        if(t > end_time && on == true)
        {
            on = false;
            start_time += cycle_length;
            end_time += cycle_length;
        }
        // If the current time is in the stimulation interval
        // then we need to activate the pacing
        else if(t >= start_time && t <= end_time)
        {
            on = true;
        }
    }
    /// True if we are applying an external stimulus
    bool on;
    /// Initial time of the current stimulus
    double start_time;
    /// Final time of the current stimulus
    double end_time;
    /// Time interval between 2 stimuli
    double cycle_length;
    /// Amplitude of the stimulus
    double amplitude;

};

/*! Computes the APD and it saves it in a vector
 *  Example: we want to measure the time interval between to V = 0.5 (APD50)
 *
 * @param on are the measurements on? This values will be changed by this function
 * @param APD_threshold percentage of the for computing the APD
 * @param time current time
 * @param APD vector to store the APDs
 * @param APD_time value where we store the last APD time
 * @param v transmembrane potential
 */
void checkAPD(bool& on, double APD_threshold, double time, std::vector<double>& APD, double& APD_time, double v)
{
		if(on)
		{
			if( v < APD_threshold)
			{
				APD_time -= time;
				APD.push_back(-APD_time);
				APD_time = 0.0;
				on = false;
			}
		}
		else
		{
		    // If the potential is greater than the given threshold
		    // we start measuring from this time
		    // Example: For APD50, if V >= 0.5 than we turn measurments on
			if( v >= APD_threshold)
			{
				on = true;
				APD_time = time;
			}
		}
}



// Start the main Program
int main(int argc, char** argv)
{
    // Print BeatIt on screen
    BeatIt::printBanner(std::cout);

    // Create a the Nas-Panfilov model using the Ionic Model factory
	std::unique_ptr<BeatIt::IonicModel> pNP( BeatIt::IonicModel::IonicModelFactory::Create("Grandi11") );
	// Check the number of variables in the model
	int numVar = pNP->numVariables();
	// Create a container to store the variables of the ionic model
	std::vector<double> variables(numVar, 0.0);

	// Read the input file
	GetPot commandLine(argc, argv);
    std::string datafile_name = commandLine.follow ( "data.beat", 2, "-i", "--input" );
    GetPot data(datafile_name);


    // Create the output file
    std::string output_file = data("output", "results.txt");
	std::ofstream output(output_file);
    // Initialize the output file
	pNP->initializeSaveData(output);
    // set the default initial conditions in the ionic model
	pNP->setup(data, "model");
	bool rv = data("model/Grandi11.resting_values", false);
	std::cout << "RV: " << rv << std::endl;
	pNP->initialize(variables);


    // We check the APD90
    bool on90 = false;
    std::vector<double> APD90;
    double APD90_time = 0.0;
    // the thershold is 0.1  because we want 90% of the action potential
    double APD90_threshold = data("threshold90", -60);
    bool on30 = false;
    std::vector<double> APD30;
    double APD30_time = 0.0;
    // the thershold is 0.1  because we want 90% of the action potential
    double APD30_threshold = data("threshold30", -30);


    // timestep
	double dt = data("dt", 0.005);
    // output at every save_iter
	int save_iter = data("savedt", 1.0) / dt;
	// Final time
	double TF = data("end", 3000.);
	// External Stimulus
	double Ist = 0;
	// Current time
	double time = 0.0;
	// save to output the initial data
	BeatIt::saveData(time, variables, output);
	// I don't remember why we use this (may be we do not use it ...)
	int iter = 0;


    // Create the pacing object and set the pacing conditions
	Pacing pacing;
	// Amplitude of the stimulus Iapp = 1
	pacing.amplitude = data("amplitude", 80.0);
	// Start to stimulate at t = 0
    pacing.start_time= data("ti", 0.0);
    // End the stimulation at t = 0.5
    pacing.end_time= pacing.start_time + data("duration", 0.5);
    // Stimulate every 15
    pacing.cycle_length= data("cl", 1000);
    pacing.on = false;


    // loop over time
    while( time <= TF )
	{
        // compute external stimulus
	    Ist = pacing.istim(time);
	    // solve one timestep
		pNP->solve(variables, Ist, dt);
		// update the current time
		time += dt;
		// update current iteration
		++iter;
		// if we want to output we can
		if( 0 == iter%save_iter ) BeatIt::saveData(time, variables, output);
		// update the APD
		checkAPD(on90, APD90_threshold, time, APD90, APD90_time, variables[0]);
		checkAPD(on30, APD30_threshold, time, APD30, APD30_time, variables[0]);
	}
    // close output file
	output.close();

	//up to the 16th digit
	const double reference_solution_norm = 25.3500000005533;

	std::cout << "APD90 size: " << APD90.size() << std::endl;
	//std::cout << std::setprecision(15) << APD90[APD90.size()-1] << std::endl;

    std::string output_apd90 = data("outputAPD90", "APD90.txt");
	std::ofstream outputAPD90(output_apd90);
	for(auto && apd90 : APD90) outputAPD90 << apd90 << std::endl;
    std::string output_apd30 = data("outputAPD30", "APD30.txt");
	std::ofstream outputAPD30(output_apd30);
	for(auto && apd30 : APD30) outputAPD30 << apd30 << std::endl;

	return BeatIt::CTest::check_test(APD90[APD90.size()-1], reference_solution_norm, 1e-12);
}

