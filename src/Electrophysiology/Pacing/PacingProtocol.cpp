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
 * \file PacingProtocol.cpp
 *
 * \class PacingProtocol
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
 * Created on: Aug 21, 2016
 *
 */

#include "Electrophysiology/Pacing/PacingProtocol.hpp"
#include "libmesh/function_base.h"
#include "libmesh/point.h"

namespace BeatIt
{


bool isPointInside(DistanceType type, double r, double x, double y ,  double z)
{
	bool isInside = false;
	switch(type)
	{
		case DistanceType::l_1:
		{
			if( r >=  std::abs(x) + std::abs(y) + std::abs(z)  ) isInside = true;
			break;
		}
		case DistanceType::l_inf:
		{
			double max = std::max(std::abs(x), std::abs(y));
			max = std::max(max, std::abs(z));
			if( r >=  max) isInside = true;
			break;
		}
		default:
		case DistanceType::l_2:
		{
			if( r*r >= x*x + y*y + z*z  ) isInside = true;
			break;
		}
        case DistanceType::l_2_2Dx:
        {
            if( r*r >= y*y + z*z  ) isInside = true;
            break;
        }
        case DistanceType::l_2_2Dy:
        {
            if( r*r >= x*x + z*z  ) isInside = true;
            break;
        }
        case DistanceType::l_2_2Dz:
        {
            if( r*r >= x*x + y*y ) isInside = true;
            break;
        }

	}
	return isInside;
}

PacingProtocol::PacingProtocol()
    : M_pacing(nullptr)
    , M_isPacingOn(false)
    , M_amplitude(10.0)
    , M_stopTime(-1.0)
    , M_startTime(0.0)
    , M_endTime(2.0)
    , M_radius(1.0)
    , M_type(DistanceType::l_2)
    , M_x0(0.0)
    , M_y0(0.0)
    , M_z0(0.0)
    , M_duration(2.0)
    , S_beats(0)
{
}

PacingProtocol::~PacingProtocol()
{
//    if(M_pacing) delete M_pacing;
}



double
PacingProtocol::eval(const Point & p,
                           const double time)
{
    double pacing = 0.0;
    bool ispInside = BeatIt::isPointInside(M_type, M_radius, p(0)-M_x0, p(1)-M_y0, p(2)-M_z0);
    if(ispInside)
    {
        if(M_isPacingOn) pacing = M_amplitude;
    }

    if(  M_stopTime > 0 && time > M_stopTime) pacing = 0.0;
    return pacing;
}


void
PacingProtocol::set_distance_type(std::string& type)
{
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
    else if(type == "l_2_2Dx")
    {
        M_type = DistanceType::l_2_2Dx;
        std::cout << "\t\t type: cylinder (x)" << std::endl;
    }
    else if(type == "l_2_2Dy")
    {
        M_type = DistanceType::l_2_2Dy;
        std::cout << "\t\t type: cylinder (y)" << std::endl;
    }
    else if(type == "l_2_2Dz")
    {
        M_type = DistanceType::l_2_2Dz;
        std::cout << "\t\t type: cylinder (z)" << std::endl;
    }
    else
    {
        M_type = DistanceType::l_2;
        std::cout << "\t\t type: sphere" << std::endl;
    }
}


} /* namespace BeatIt */
