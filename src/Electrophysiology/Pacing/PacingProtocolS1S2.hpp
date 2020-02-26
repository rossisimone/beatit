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
 * PacingProtocolS1S2.hpp
 *
 *  Created on: Nov 2, 2016
 *      Author: srossi
 */

#ifndef SRC_ELECTROPHYSIOLOGY_PACING_PacingProtocolS1S2_HPP_
#define SRC_ELECTROPHYSIOLOGY_PACING_PacingProtocolS1S2_HPP_

#include "Electrophysiology/Pacing/PacingProtocol.hpp"
#include "libmesh/function_base.h"

namespace BeatIt {


class PacingProtocolS1S2 : public virtual PacingProtocol,
											public virtual libMesh::FunctionBase<libMesh::Number>
{
public:
	 typedef libMesh::FunctionBase<libMesh::Number> super;
	PacingProtocolS1S2();
	virtual ~PacingProtocolS1S2();
	void setup(const GetPot& data, std::string section = "monodomain/pacing");

	void update(double time);
	    /**
     * Returns a new copy of the function.  The new copy should be as
     * ``deep'' as necessary to allow independent destruction and
     * simultaneous evaluations of the copies in different threads.
     */
    std::unique_ptr<libMesh::FunctionBase<double> > clone () const;


    void operator() (const  Point & p,
                    const double time,
                    libMesh::DenseVector<double> & output);
    double component(unsigned int i,
                    const Point & p,
                    double time = 0.0);
    double operator() (const Point & p,
                              const double time = 0.);


    double  M_S1cycleLength;
    double  M_S2cycleLength;
    int     M_numS1Stimuli;
    int     M_numS2Stimuli;
    int     M_S1S2Counter;
    double  M_minCycleLength;
    double M_S1Decrement;
    double M_S2Decrement;
    unsigned int M_decrementBeats;
};

PacingProtocol* createPacingProtocolS1S2();

namespace
{
	static bool register_PacingProtocolS1S2 = BeatIt::PacingProtocol::PacingProtocolFactory::Register("S1S2", &createPacingProtocolS1S2);
}


} /* namespace BeatIt */

#endif /* SRC_ELECTROPHYSIOLOGY_PACING_PacingProtocolS1S2_HPP_ */
