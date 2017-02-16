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
 * PacingProtocolS1.cpp
 *
 *  Created on: Nov 2, 2016
 *      Author: srossi
 */

#include "Electrophysiology/Pacing/PacingProtocolS1.hpp"
#include "libmesh/getpot.h"
#include "libmesh/point.h"

namespace BeatIt {


PacingProtocol* createPacingProtocolS1()
{
	return new PacingProtocolS1;
}
PacingProtocolS1::PacingProtocolS1()
    :  M_isPacingOn(false)
    , M_startTime(0.0)
    , M_endTime(2.0)
    , M_cycleLength(60.0)
    , M_minCycleLength(10.0)
    , M_amplitude(10.)
    , M_type(DistanceType::l_2)
    , M_radius(1.0)
    , M_decrement(0.0)
    , M_decrementBeats(5)
	, M_x0(0.0)
	, M_y0(0.0)
	, M_z0(0.0)
{
	_initialized = true;
    _is_time_dependent = true;
}

PacingProtocolS1::~PacingProtocolS1()
{
}

void
PacingProtocolS1::setup(const GetPot& data, std::string section)
{
    std::cout << "* PacingProtocolS1: reading  from input file" << std::endl;
    M_cycleLength  = data( section+"/cycle_length", 60.);
    M_amplitude  = data( section+"/amplitude", 10.);
    M_startTime  = data( section+"/start_time", 0.0);
    M_endTime = M_startTime + data( section+"/duration", 1.0);
    M_radius = data( section+"/radius", 0.15);
    M_decrement = data( section+"/decrement", 0.0);
    M_decrementBeats = data( section+"/decrement_beats", 5);
    M_minCycleLength = data( section+"/cycle_length_min", 10.0);
    M_stopTime = data( section+"/stop_time", -1.0);
    M_x0 = data( section+"/x0", 0.0);
    M_y0 = data( section+"/y0", 0.0);
    M_z0 = data( section+"/z0", 0.0);
    std::cout << "\t\t start_time: " << M_startTime << std::endl;
    std::cout << "\t\t end_time: " << M_endTime << std::endl;
    std::cout << "\t\t cycle_length: " << M_cycleLength << std::endl;
    std::cout << "\t\t radius: " << M_radius << std::endl;
    std::cout << "\t\t center: " << M_x0 << ", " << M_y0 << ", " << M_z0 << std::endl;
    std::cout << "\t\t amplitude: " << M_amplitude << std::endl;

	std::string type = data(section+"/distance", "l_2");
    if(type == "l_1")
	{
    	M_type = DistanceType::l_1;
    	std::cout << "\t\t type: diamond" << std::endl;
	}
    else if(type == "l_inf")
	{
    	M_type = DistanceType::l_inf;
    	std::cout << "\t\t type: box" << std::endl;
	}
    else
    {
    	    	std::cout << "\t\t type: sphere" << std::endl;
    }
	M_pacing = this;
	std::cout << "\t\t this: " << this << ", pacing: " << M_pacing << std::endl;

}



	/**
 * Returns a new copy of the function.  The new copy should be as
 * ``deep'' as necessary to allow independent destruction and
 * simultaneous evaluations of the copies in different threads.
 */
libMesh::UniquePtr<libMesh::FunctionBase<double> >
PacingProtocolS1::clone () const
{
	PacingProtocolS1* pcopy = new PacingProtocolS1;
	pcopy->_master = _master;
    pcopy->_is_time_dependent = true;
    pcopy->_initialized = true;
    pcopy-> M_isPacingOn = M_isPacingOn;
    pcopy-> M_startTime = M_startTime;
    pcopy-> M_endTime = M_endTime;
    pcopy-> M_cycleLength = M_cycleLength;
    pcopy-> M_amplitude = M_amplitude;
    pcopy-> M_type = M_type;
    pcopy-> M_radius = M_radius;
    pcopy-> M_decrementBeats = M_decrementBeats;
    pcopy-> M_minCycleLength = M_minCycleLength;
    pcopy-> M_decrement = M_decrement;
    pcopy-> M_x0 = M_x0;
    pcopy-> M_y0 = M_y0;
    pcopy-> M_z0 = M_z0;
    pcopy-> M_stopTime = M_stopTime;

	return libMesh::UniquePtr<libMesh::FunctionBase<double> >(pcopy);
}


double
PacingProtocolS1::operator() (const Point & p,
						   const double time)
{
	double pacing = 0.0;
	if(time < M_startTime) M_isPacingOn = false;
	else if(time >= M_startTime && time <= M_endTime) M_isPacingOn = true;
	else // time > M_endTime
	{
	    S_beats++;
	    if(S_beats%M_decrementBeats == 0  && M_cycleLength > M_minCycleLength)
        {
	        M_cycleLength -= M_decrement;
        }
		M_startTime += M_cycleLength;
		M_endTime += M_cycleLength;
    	M_isPacingOn = false;
	}

	bool isPointInside = BeatIt::isPointInside(M_type, M_radius, p(0)-M_x0, p(1)-M_y0, p(2)-M_z0);
	if(isPointInside)
	{
		if(M_isPacingOn) pacing = M_amplitude;
	}

	if(  M_stopTime > 0 && time > M_stopTime) pacing = 0.0;
	return pacing;
}
void
PacingProtocolS1::operator() (const Point & p,
													 const double time,
													 libMesh::DenseVector<double> & output)
{
	output(0) =   operator()(p, time);
}
double
PacingProtocolS1::component(unsigned int i,
															 const Point & p,
															 double time)
{
	return operator()(p, time);
}

} /* namespace BeatIt */
