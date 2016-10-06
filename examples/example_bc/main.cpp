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
 * main.cpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */

#include "BoundaryConditions/BCHandler.hpp"
#include "libmesh/getpot.h"
#include "BoundaryConditions/BCData.hpp"
int main (int argc, char ** argv)
{
	GetPot data("data.pot");

	BeatIt::BCHandler bc1;
	bc1.readBC(data,"physic1");
	BeatIt::BCHandler bc2;
	bc2.readBC(data, "physic2");
	bc1.showMe();
	bc2.showMe();

	double x = 1.0;
	double y = 2.0;
	double z = 3.0;
	double t = 0.0;

    // bc1 function: 1.0
	double bc1_val_0 = (bc1.get_bc(1) ->get_function())(t,x,y,z,0);
    std::cout << "bc1_val_0 = " << bc1_val_0 << std::endl;
    double bc1_val_1 = (bc1.get_bc(1) ->get_function())(t,x,y,z,1);
    std::cout << "bc1_val_1 = " << bc1_val_1 << std::endl;
    double bc1_val_2 = (bc1.get_bc(1) ->get_function())(t,x,y,z,2);
    std::cout << "bc1_val_2 = " << bc1_val_2 << std::endl;
    double bc1_val      = (bc1.get_bc(2) ->get_function())(t,x,y,z,0);
    std::cout << "bc1_val = " << bc1_val << std::endl;
    double bc2_val      = (bc2.get_bc(3) ->get_function())(t,x,y,z,0);
    std::cout << "bc2_val = " << bc2_val << std::endl;

	return 0;
}

