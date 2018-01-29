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
 * BCData.cpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */

#include "BoundaryConditions/BCData.hpp"
#include "libmesh/getpot.h"
#include "Util/IO/io.hpp"

namespace BeatIt {


BCData::ModeMap  BCData::S_modeMap =
{
		{"Full", BCMode::Full},
		{"Component", BCMode::Component},
		{"Normal", BCMode::Normal},
		{"Tangential", BCMode::Tangential}
};

BCData::ComponentMap  BCData::S_componentMap =
{
		{"X", BCComponent::X},
		{"Y", BCComponent::Y},
		{"Z", BCComponent::Z},
		{"All", BCComponent::All}
};

BCData::TypeMap  BCData::S_typeMap =
{
		{"Dirichlet", BCType::Dirichlet},
		{"Neumann", BCType::Neumann},
		{"Robin", BCType::Robin},
        {"NitscheSymmetric", BCType::NitscheSymmetric},
        {"NitscheUnsymmetric", BCType::NitscheUnsymmetric},
        {"Penalty", BCType::Penalty}
};


BCData::BCData()
{
	// TODO Auto-generated constructor stub

}

BCData::~BCData() {
	// TODO Auto-generated destructor stub
}

void
BCData::setup(const GetPot& data, const std::string& section )
{
	std::string flags =  data(section+"/flag", "-1");
	std::vector<int> flag_vec;
	BeatIt::readList(flags,flag_vec);

//	int flag = data(section+"/flag", -1);
	std::cout << "Section: " << section << std::endl;
	for(auto && fl : flag_vec)
	{
		if(fl >= 0 ) M_flag.push_back(  static_cast<unsigned int>(fl) );
		else
		{
			std::cout << "- Error -- BCData:: cannot read  the boundary conditions from " << section << ". BCData FLAG Error!" << std::endl;
			throw std::runtime_error("BCData FLAG Error!");
		}
	}

	std::string comp = data(section+"/component", "");
	auto itc = S_componentMap.find(comp);
	if(itc != S_componentMap.end() ) M_component = itc->second;
	else
	{
		std::cout << "- Error -- BCData:: cannot read  the boundary conditions from " << section << ". BCData COMPONENT Error!" << std::endl;
		throw std::runtime_error("BCData COMPONENT Error!");
	}

	std::string type = data(section+"/type", "");
	auto itt = S_typeMap.find(type);
	if(itt != S_typeMap.end() ) M_type = itt->second;
	else
	{
		std::cout << "- Error -- BCData:: cannot read  the boundary conditions from " << section << ". BCData TYPE Error!" << std::endl;
		throw std::runtime_error("BCData TYPE Error!");
	}

	std::string mode = data(section+"/mode", "");
	auto itm = S_modeMap.find(mode);
	if(itm != S_modeMap.end() ) M_mode = itm->second;
	else
	{
		std::cout << "- Error -- BCData:: cannot read  the boundary conditions from " << section << ". BCData MODE Error!" << std::endl;
		throw std::runtime_error("BCData MODE Error!");
	}

	std::string function = data(section+"/function", "");
	if("fe_function" == function)
	{
	    M_using_fe_function = true;
	}
	else if( "" != function)
	{
		M_function.read(function);
	}
	else
	{
		std::cout << "- Error -- BCData:: cannot read  the boundary conditions from " << section << ". BCData FUNCTION Error!" << std::endl;
		throw std::runtime_error("BCData FUNCTION Error!");
	}

}


void
BCData::showMe( std::ostream& ofstream  )
{
	ofstream << " \t BCData: " << std::endl;
	std::cout << "\t\t Num flags: " << M_flag.size() << std::endl;
	ofstream << " \t\t FLAGs: ";
	for (auto && fl : M_flag ) std::cout << fl  << ", ";
	std::cout << std::endl;
	ofstream << " \t\t TYPE: ";
	switch(M_type)
	{
		case BCType::Dirichlet:
		{
			ofstream << " Dirichlet" << std::endl;
			break;
		}
		case BCType::Neumann:
		{
			ofstream << " Neumann" << std::endl;
			break;
		}
		case BCType::Robin:
		{
			ofstream << " Robin" << std::endl;
			break;
		}
        case BCType::NitscheSymmetric:
        {
            ofstream << " NitscheSymmetric" << std::endl;
            break;
        }
        case BCType::NitscheUnsymmetric:
        {
            ofstream << " NitscheUnsymmetric" << std::endl;
            break;
        }
        case BCType::Penalty:
        {
            ofstream << " Penalty" << std::endl;
            break;
        }
		default:
		{
			throw std::runtime_error("BCData: Showme - Wrong Type");
			break;
		}
	}
	ofstream << " \t\t MODE: ";
	switch(M_mode)
	{
		case BCMode::Full:
		{
			ofstream << " Full" << std::endl;
			break;
		}
		case BCMode::Component:
		{
			ofstream << " Component" << std::endl;
			break;
		}
		case BCMode::Normal:
		{
			ofstream << " Normal" << std::endl;
			break;
		}
		case BCMode::Tangential:
		{
			ofstream << " Tangential" << std::endl;
			break;
		}
		default:
		{
			throw std::runtime_error("BCData: Showme - Wrong Mode");
			break;
		}
	}


	ofstream << " \t\t COMPONENT: ";
	switch(M_component)
	{
		case BCComponent::X:
		{
			ofstream << " X" << std::endl;
			break;
		}
		case BCComponent::Y:
		{
			ofstream << " Y" << std::endl;
			break;
		}
		case BCComponent::Z:
		{
			ofstream << " Z" << std::endl;
			break;
		}
		case BCComponent::All:
		{
			ofstream << " All" << std::endl;
			break;
		}
		default:
		{
			throw std::runtime_error("BCData: Showme - Wrong Component");
			break;
		}
	}

	if(M_using_fe_function)
	{
	    ofstream << " \t\t Using FE function \n";
	}
	else
	{
        ofstream << " \t\t FUNCTION: \n";
        M_function.showMe(ofstream);
	}
}

double
BCData::fe_function_component(double t, double x, double y, double z, int component)
{
    if(M_fe_function)
    {
        return M_fe_function->component(component, libMesh::Point(x,y,z), t);
    }
    else
    {
        throw std::runtime_error("Calling fe_function_component, the fe_function was not set!");
        return 0.0;
    }
}


} /* namespace BeatIt */
