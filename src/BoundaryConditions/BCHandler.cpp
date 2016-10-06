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
 * BCHandler.cpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */

#include "BoundaryConditions/BCHandler.hpp"
#include "BoundaryConditions/BCData.hpp"
#include "libmesh/getpot.h"
#include "Util/IO/io.hpp"
//#include <iostream>
namespace BeatIt
{

BCHandler::BCHandler()
{

}

BCHandler::~BCHandler()
{
}

void
BCHandler::readBC(const GetPot& data,  const std::string& section)
{
	std::string BCList = data(section+"/BC/list", "");
	std::vector<std::string> BCListVector;
	bool read = BeatIt::readList(BCList, BCListVector);
	if(read)
	{
		for(auto&& b : BCListVector)
		{
			std::shared_ptr<BCData> bc_ptr( new BCData() );
		    std::string path = section+"/BC/"+b;
			bc_ptr->setup(data, path);
			M_bcNames.push_back(path);
			M_bcs.push_back(bc_ptr);
			auto size = bc_ptr->size();
			if(size <= 0)
			{
				std::cout << "Wrong size in readBC" << std::endl;
				throw std::runtime_error("Error in readBC");
			}
			for(int i = 0; i < size; ++i) M_bcMap[bc_ptr->get_flag(i)] = bc_ptr;
		}
	}

}

void
BCHandler::showMe(std::ostream& ofstream )
{
	ofstream << "\t\t\t BCHANDLER show me: " << std::endl;
	int size = M_bcs.size();
	for (int i = 0; i < size; ++i)
	{
		ofstream << "\t "<< M_bcNames[i] << std::endl;
		M_bcs[i]->showMe(ofstream);
	}
	std::cout << std::endl;
}

std::shared_ptr<BCData>
BCHandler::get_bc(unsigned int ID)
{
	auto it = M_bcMap.find(ID);
	if( it != M_bcMap.end() ) return it->second;
	else
	{
		std::shared_ptr<BCData> ptr;
		return ptr;
	}
}



} // BeatIt
