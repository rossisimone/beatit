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

namespace BeatIt
{

int PacingProtocol::S_beats = 0;

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
	}
	return isInside;
}

PacingProtocol::PacingProtocol()
    : M_pacing(nullptr)
{
}

PacingProtocol::~PacingProtocol()
{
//    if(M_pacing) delete M_pacing;
}



} /* namespace BeatIt */
