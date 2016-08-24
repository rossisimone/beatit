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
 * \file TimeData.cpp
 *
 * \class TimeData
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
 * Created on: Aug 11, 2016
 *
 */

#include "Util/TimeData.hpp"
#include <iostream>
namespace BeatIt
{

TimeData::TimeData()
    : M_dt(1.0)
    , M_startTime(0.0)
    , M_endTime(1.0)
    , M_maxIter(1000000)
    , M_iter(0)
{
}


void
TimeData::setup(const GetPot& data, const std::string& section)
{
    M_dt = data(section+"time/dt", 1.0);
    M_startTime = data(section+"time/init_time", 0.0);
    M_time = M_startTime;
    M_endTime = data(section+"time/final_time", 1.0);
    M_maxIter = data(section+"time/max_iter", 99999999);
    M_saveIter = data(section+"time/save_iter",1);
}

void
TimeData::print()
{
    std::cout << "* TimeData parameters: " << std::endl;
    std::cout << "\tdt        = " << M_dt << std::endl;
    std::cout << "\ttime      = " << M_time << std::endl;
    std::cout << "\tstart     = " << M_startTime << std::endl;
    std::cout << "\tend       = " << M_endTime << std::endl;
    std::cout << "\titer      = " << M_iter << std::endl;
    std::cout << "\tmax_iter  = " << M_maxIter << std::endl;
    std::cout << "\tsave_iter = " << M_saveIter << std::endl;
}


} /* namespace BeatIt */
