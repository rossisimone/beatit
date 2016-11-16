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
 * PacingProtocolS1.hpp
 *
 *  Created on: Nov 2, 2016
 *      Author: srossi
 */

#ifndef SRC_ELECTROPHYSIOLOGY_PACING_PACINGPROTOCOLS1_HPP_
#define SRC_ELECTROPHYSIOLOGY_PACING_PACINGPROTOCOLS1_HPP_

#include "Electrophysiology/Pacing/PacingProtocol.hpp"
#include "libmesh/function_base.h"

namespace BeatIt {


class PacingProtocolS1 : public virtual PacingProtocol,
											public virtual libMesh::FunctionBase<libMesh::Number>
{
public:
	 typedef libMesh::FunctionBase<libMesh::Number> super;
	     typedef libMesh::Point Point;
	PacingProtocolS1();
	virtual ~PacingProtocolS1();
	void setup(const GetPot& data, std::string section = "monodomain/pacing");


	    /**
     * Returns a new copy of the function.  The new copy should be as
     * ``deep'' as necessary to allow independent destruction and
     * simultaneous evaluations of the copies in different threads.
     */
    libMesh::UniquePtr<libMesh::FunctionBase<double> > clone () const;


    double operator() (const Point & p,
                               const double time = 0.);
    void operator() (const Point & p,
                     const double time,
                     libMesh::DenseVector<double> & output);
    double component(unsigned int i,
                     const Point & p,
                     double time = 0.0);

    bool M_isPacingOn;
    double  M_startTime;
    double  M_endTime;
    double  M_cycleLength;
    double M_minCycleLength;
    double M_amplitude;
    double M_radius;
    DistanceType M_type;
    double M_decrement;
    unsigned int M_decrementBeats;
    double M_x0;
    double M_y0;
    double M_z0;
};

PacingProtocol* createPacingProtocolS1();

namespace
{
	static bool register_PacingProtocolS1 = BeatIt::PacingProtocol::PacingProtocolFactory::Register("S1", &createPacingProtocolS1);
}


} /* namespace BeatIt */

#endif /* SRC_ELECTROPHYSIOLOGY_PACING_PACINGPROTOCOLS1_HPP_ */
