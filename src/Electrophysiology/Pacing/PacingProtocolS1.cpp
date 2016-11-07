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
    , M_cycleLength(10.0)
    , M_amplitude(10.)
    , M_type(DistanceType::l_2)
    , M_radius(1.0)
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
    M_cycleLength  = data( section+"/cycle_length", 10.);
    M_amplitude  = data( section+"/amplitude", 10.);
    M_startTime  = data( section+"/start_time", 0.0);
    M_endTime = M_startTime + data( section+"/duration", 1.0);;
    M_radius = data( section+"/radius", 0.15);

    std::cout << "\t\t start_time: " << M_startTime << std::endl;
    std::cout << "\t\t end_time: " << M_endTime << std::endl;
    std::cout << "\t\t cycle_length: " << M_cycleLength << std::endl;
    std::cout << "\t\t radius: " << M_radius << std::endl;
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
		M_startTime += M_cycleLength;
		M_endTime += M_cycleLength;
    	M_isPacingOn = false;
	}

	bool isPointInside = BeatIt::isPointInside(M_type, M_radius, p(0), p(1), p(2));
	if(isPointInside)
	{
		if(M_isPacingOn) pacing = M_amplitude;
	}
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
