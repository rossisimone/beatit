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
 * \file MonodomainUtil.cpp
 *
 *
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
 * Created on: Aug 15, 2016
 *
 */



#include "libmesh/parallel.h"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"
#include "libmesh/numeric_vector.h"
#include "libmesh/explicit_system.h"


namespace BeatIt
{

namespace Electrophysiology
{


bool isFullyActivated(libMesh::System& system, double initValue )
{

	int local = 1;
    auto first_local_index = system.solution->first_local_index();
    auto last_local_index = system.solution->last_local_index();
    for(auto i = first_local_index; i < last_local_index; i++)
    {
        double val_at =  (*system.solution)(i);
    	if(initValue ==  val_at )
    	{
    		local = 0;
    		break;
    	}
    }

    system.get_mesh().comm().min<int>(local);
    if( 0 == local ) return false;
    else return true;
}



void
updateActivationTimes( const libMesh::NumericVector<libMesh::Number>& potential,
														 libMesh::System& system,
														 double time,
														 double threshold,
														 double initiValue)
{
    auto first_local_index = system.solution->first_local_index();
    auto last_local_index = system.solution->last_local_index();

    double val_at = -1.0;
    double val_v  = 0.0;

    for(auto i = first_local_index; i < last_local_index; i++)
    {
        val_at =  (*system.solution)(i);
        val_v  =  potential(i);
        if( val_at == -1.0 && val_v >= threshold )
        {
        	system.solution->set(i, time);
        }
    }
    system.solution->close();
}




} // Electrophysiology

} // BeatIt


