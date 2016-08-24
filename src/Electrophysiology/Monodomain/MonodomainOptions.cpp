/*
 * MonodomainOptions.cpp
 *
 *  Created on: Aug 15, 2016
 *      Author: srossi
 */



#include "Electrophysiology/Monodomain/MonodomainOptions.hpp"

namespace BeatIt
{


namespace Electrophysiology
{

MassMatrixType selectMassMatrixType(std::string type)
{
	std::map<std::string, MassMatrixType> map_mass;
	map_mass["consistent"] = MassMatrixType::CONSISTENT;
	map_mass["lumped"]      = MassMatrixType::LUMPED;
	map_mass["highorder"]  = MassMatrixType::HIGHORDER;
	auto it = map_mass.find(type);
	if( it != map_mass.end() ) return it->second;
	else return MassMatrixType::CONSISTENT;
}


AssemblyType    selectAssemblyType(std::string type)
{
	std::map<std::string, AssemblyType> map_assembly;
	map_assembly["matrix"] = AssemblyType::MATRIX;
	map_assembly["full"]      = AssemblyType::FULL;
	auto it = map_assembly.find(type);
	if( it != map_assembly.end() ) return it->second;
	else return AssemblyType::MATRIX;
}


} // Electrophysiology

} // BeatIt
