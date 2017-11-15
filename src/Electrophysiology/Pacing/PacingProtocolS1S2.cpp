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

/*
 * PacingProtocolS1S2.cpp
 *
 *  Created on: Nov 2, 2016
 *      Author: srossi
 */

#include "Electrophysiology/Pacing/PacingProtocolS1S2.hpp"
#include "libmesh/getpot.h"

namespace BeatIt {


PacingProtocol* createPacingProtocolS1S2()
{
	return new PacingProtocolS1S2;
}
PacingProtocolS1S2::PacingProtocolS1S2()
    : M_S1cycleLength(60.0)
    , M_S2cycleLength(60.0)
    , M_numS1Stimuli(1)
    , M_numS2Stimuli(1)
    , M_minCycleLength(10.0)
    , M_S1Decrement(0.0)
    , M_S2Decrement(0.0)
    , M_decrementBeats(5)
    , M_S1S2Counter(0)
{
	_initialized = true;
    _is_time_dependent = true;
}

PacingProtocolS1S2::~PacingProtocolS1S2()
{
}

void
PacingProtocolS1S2::setup(const GetPot& data, std::string section)
{
    std::cout << "* PacingProtocolS1S2: reading  from input file" << std::endl;
    M_S1cycleLength  = data( section+"/s1_cycle_length", 60.);
    M_S2cycleLength  = data( section+"/s2_cycle_length", 60.);
    M_amplitude  = data( section+"/amplitude", 10.);
    M_startTime  = data( section+"/start_time", 0.0);
    M_duration = data( section+"/duration", 1.0);
    M_endTime = M_startTime + M_duration;
    M_radius = data( section+"/radius", 0.15);
    M_S1Decrement = data( section+"/s1_decrement", 0.0);
    M_S2Decrement = data( section+"/s2_decrement", 0.0);
    M_decrementBeats = data( section+"/decrement_beats", 5);
    M_minCycleLength = data( section+"/cycle_length_min", 10.0);
    M_stopTime = data( section+"/stop_time", -1.0);
    M_numS1Stimuli = data( section+"/num_s1", 1);
    M_numS2Stimuli = data( section+"/num_s2", 1);
    M_boundaryID = data( section+"/boundary", -1);

    M_x0 = data( section+"/x0", 0.0);
    M_y0 = data( section+"/y0", 0.0);
    M_z0 = data( section+"/z0", 0.0);
    std::cout << "\t\t start_time: " << M_startTime << std::endl;
    std::cout << "\t\t end_time: " << M_endTime << std::endl;
    std::cout << "\t\t s1_cycle_length: " << M_S1cycleLength  << ", num_s1: " << M_numS1Stimuli << std::endl;
    std::cout << "\t\t s2_cycle_length: " << M_S2cycleLength  << ", num_s2: " << M_numS2Stimuli << std::endl;
    std::cout << "\t\t radius: " << M_radius << std::endl;
    std::cout << "\t\t center: " << M_x0 << ", " << M_y0 << ", " << M_z0 << std::endl;
    std::cout << "\t\t amplitude: " << M_amplitude << std::endl;
    std::cout << "\t\t boundary: " << M_boundaryID << std::endl;
    std::string type = data(section+"/distance", "l_2");
    PacingProtocol::set_distance_type(type);

	M_pacing = this;
	std::cout << "\t\t this: " << this << ", pacing: " << M_pacing << std::endl;

}



	/**
 * Returns a new copy of the function.  The new copy should be as
 * ``deep'' as necessary to allow independent destruction and
 * simultaneous evaluations of the copies in different threads.
 */
libMesh::UniquePtr<libMesh::FunctionBase<double> >
PacingProtocolS1S2::clone () const
{
	PacingProtocolS1S2* pcopy = new PacingProtocolS1S2;
	pcopy->_master = _master;
    pcopy->_is_time_dependent = true;
    pcopy->_initialized = true;
    pcopy-> M_isPacingOn = M_isPacingOn;
    pcopy-> M_startTime = M_startTime;
    pcopy-> M_endTime = M_endTime;
    pcopy-> M_S1cycleLength = M_S1cycleLength;
    pcopy-> M_S2cycleLength = M_S2cycleLength;
    pcopy-> M_numS1Stimuli = M_numS1Stimuli;
    pcopy-> M_numS2Stimuli = M_numS2Stimuli;
    pcopy-> M_amplitude = M_amplitude;
    pcopy-> M_type = M_type;
    pcopy-> M_radius = M_radius;
    pcopy-> M_decrementBeats = M_decrementBeats;
    pcopy-> M_minCycleLength = M_minCycleLength;
    pcopy-> M_S1Decrement = M_S1Decrement;
    pcopy-> M_S2Decrement = M_S2Decrement;
    pcopy-> M_S1S2Counter = M_S1S2Counter;
    pcopy-> M_x0 = M_x0;
    pcopy-> M_y0 = M_y0;
    pcopy-> M_z0 = M_z0;
    pcopy-> M_stopTime = M_stopTime;

	return libMesh::UniquePtr<libMesh::FunctionBase<double> >(pcopy);
}


void
PacingProtocolS1S2::update(double time)
{
	if(time < M_startTime) M_isPacingOn = false;
	else if(time >= M_startTime && time <= M_endTime)
    {
	    M_isPacingOn = true;
    }
	else // time > M_endTime
	{
	    M_S1S2Counter++;
	    S_beats++;
	    if(S_beats%M_decrementBeats == 0  && M_S2cycleLength > M_minCycleLength)
        {
	        M_S1cycleLength -= M_S1Decrement;
            M_S2cycleLength -= M_S2Decrement;
        }
	    if(M_S1S2Counter >= M_numS1Stimuli)
	    {
            M_startTime += M_S2cycleLength;
            M_endTime += M_S2cycleLength;
	    }
	    else
	    {
            M_startTime += M_S1cycleLength;
            M_endTime += M_S1cycleLength;
	    }
	    if(M_S1S2Counter >= M_numS1Stimuli + M_numS1Stimuli) M_S1S2Counter = 0;
    	M_isPacingOn = false;
	}
}

void
PacingProtocolS1S2::operator() (const Point & p,
                             const double time,
                             libMesh::DenseVector<double> & output)
{
    output(0) =   operator()(p, time);
}

double
PacingProtocolS1S2::component(unsigned int i,
                         const Point & p,
                         double time)
{
    return operator()(p, time);
}

double
PacingProtocolS1S2::operator() (const Point & p,
                          const double time)
{
    return PacingProtocol::eval(p, time);
}



} /* namespace BeatIt */
