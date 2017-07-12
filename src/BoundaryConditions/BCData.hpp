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
 * BCData.hpp
 *
 *  Created on: Sep 14, 2016
 *      Author: srossi
 */

#ifndef SRC_BOUNDARYCONDITIONS_BCDATA_HPP_
#define SRC_BOUNDARYCONDITIONS_BCDATA_HPP_

#include "Util/SpiritFunction.hpp"
#include <fstream>

class GetPot;

namespace BeatIt {


enum class BCMode { Full, Component, Normal, Tangential };
enum class BCComponent{ X, Y,  Z, All };
enum class BCType   { Dirichlet, Neumann, Robin, NitscheSymmetric, NitscheUnsymmetric, Penalty };

class BCData {
public:
	BCData();
	virtual ~BCData();
	void setup(const GetPot& data, const std::string& section = "" );
	void showMe(std::ostream& ofstream = std::cout );

	unsigned int size() const
	{
		return M_flag.size();
	}

	unsigned int get_flag(int index = 0) const
	{
		return M_flag[index];
	}

	BCType get_type() const
	{
		return M_type;
	}
	BCMode get_mode() const
	{
		return M_mode;
	}
	 BCComponent get_component() const
	{
		return M_component;
	}

	const SpiritFunction&  get_function() const
	{
		return M_function;
	}

protected:
	std::vector<unsigned int>         M_flag;
	SpiritFunction          M_function;
	BCComponent        M_component;
    BCMode                   M_mode;
    BCType                    M_type;

    typedef std::map<std::string, BCMode> ModeMap;
    typedef std::map<std::string, BCComponent> ComponentMap;
    typedef std::map<std::string, BCType> TypeMap;
    static ModeMap S_modeMap;
    static ComponentMap S_componentMap;
    static TypeMap S_typeMap;

};

} /* namespace BeatIt */

#endif /* SRC_BOUNDARYCONDITIONS_BCDATA_HPP_ */
