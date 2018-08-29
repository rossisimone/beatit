/*
 * SetNumericVectorValues.hpp
 *
 *  Created on: May 14, 2018
 *      Author: srossi
 */

#ifndef SRC_UTIL_SETNUMERICVECTORVALUES_HPP_
#define SRC_UTIL_SETNUMERICVECTORVALUES_HPP_


namespace libMesh
{
 class EquationSystems;
 class System;
}

namespace BeatIt
{

namespace SetValues
{

    void set_scalar_solution_on_boundary(libMesh::EquationSystems &es, libMesh::System & system, unsigned int boundID, double value, int subdomain);


} // SetValues

}// BeatIt



#endif /* SRC_UTIL_SETNUMERICVECTORVALUES_HPP_ */
