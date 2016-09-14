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
 * BCHandler.hpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */

#ifndef SRC_BOUNDARYCONDITIONS_BCHANDLER_HPP_
#define SRC_BOUNDARYCONDITIONS_BCHANDLER_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <memory>

class GetPot;

namespace BeatIt
{

class BCData;

class BCHandler {
public:
	BCHandler();
	virtual ~BCHandler();

	void readBC(GetPot& data,  const std::string& section = "");

	void showMe(std::ostream& ofstream = std::cout );

	std::vector<std::shared_ptr<BCData> > M_bcs;
	std::vector<std::string> M_bcNames;
};


} // Beat It

#endif /* SRC_BOUNDARYCONDITIONS_BCHANDLER_HPP_ */
