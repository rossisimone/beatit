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
 * \file PacingProtocol.hpp
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

#ifndef SRC_ELECTROPHYSIOLOGY_PACINGPROTOCOL_HPP_
#define SRC_ELECTROPHYSIOLOGY_PACINGPROTOCOL_HPP_

#include "libmesh/function_base.h"
#include <string>
#include <memory>
#include <cmath>
#include "Util/Factory.hpp"

class GetPot;

namespace BeatIt
{

enum class DistanceType { l_1,  l_2,  l_inf };

bool isPointInside(DistanceType type, double r, double x, double y = 0.0,  double z = 0.0);



/*!
 *
 */
class PacingProtocol
{
public:

		/// Create a factory
	typedef Factory<PacingProtocol, std::string>     PacingProtocolFactory;
    typedef libMesh::Point Point;

    PacingProtocol();
    virtual ~PacingProtocol();

    libMesh::FunctionBase<double>& pacing()
    {
        return *M_pacing;
    }

    virtual void setup(const GetPot& data, std::string section = "monodomain/pacing") = 0;
    virtual void update(double time) = 0;

    libMesh::FunctionBase<double>* M_pacing;
    double eval  (const Point & p,
                   const double time = 0.);

    static int S_beats;
    bool M_isPacingOn;
    double M_amplitude;
    double M_stopTime;
    double  M_startTime;
    double  M_endTime;
    double M_radius;
    DistanceType M_type;
    double M_x0;
    double M_y0;
    double M_z0;
    double M_duration;
};

} /* namespace BeatIt */

#endif /* SRC_ELECTROPHYSIOLOGY_PACINGPROTOCOL_HPP_ */
