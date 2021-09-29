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
 * \file TimeData.hpp
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

#ifndef SRC_UTIL_TIMEDATA_HPP_
#define SRC_UTIL_TIMEDATA_HPP_

#include "libmesh/getpot.h"
namespace BeatIt
{

/*!
 *
 */
class TimeData
{
public:
    TimeData();
    void setup(const GetPot& data, const std::string& section = "");
    void advance() { ++M_iter; M_time += M_dt; }
    void print();
    void reset();
    double M_dt;
    double M_startTime;
    double M_time;
    double M_endTime;
    int    M_maxIter;
    int    M_iter;
    int     M_saveIter;
    std::string M_section;

};

} /* namespace BeatIt */

#endif /* SRC_UTIL_TIMEDATA_HPP_ */
