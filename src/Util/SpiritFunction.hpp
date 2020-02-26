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
 * \file SpiritFunction.hpp
 *
 * \class SpiritFunction
 *
 * \brief This class provides a simple factory implementation
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
 * Created on: Aug 19, 2016
 *
 */

#ifndef SRC_UTIL_SPIRITFUNCTION_HPP_
#define SRC_UTIL_SPIRITFUNCTION_HPP_

#include "Util/Grammar.hpp"
#include "libmesh/function_base.h"

namespace BeatIt
{

/*!
 *
 */

class SpiritFunction : public virtual libMesh::FunctionBase<libMesh::Number>
{
public:
    typedef libMesh::FunctionBase<libMesh::Number> super;
    typedef libMesh::Point Point;
    typedef std::string::const_iterator iterator;
    typedef Grammar<iterator> grammar_Type;
    SpiritFunction();
    ~SpiritFunction();

//    /**
//     * Prepares a context object for use.
//     */

//    virtual void  p();

//    virtual void init_context (const FEMContext & c) libmesh_override
//    {
//        for (unsigned int v=0; v != _n_vars; ++v)
//           {
//            libMesh::FEBase * elem_fe;
//             c.get_element_fe(v, elem_fe);
//             if (_n_requested_vars)
//               elem_fe->get_phi();
//             if (_n_requested_grad_components)
//               elem_fe->get_dphi();
//       #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
//             if (_n_requested_hess_components)
//               elem_fe->get_d2phi();
//       #endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
//           }
//
//         if (_requested_normals)
//           {
//             libMesh::FEBase * side_fe;
//             c.get_side_fe(0, side_fe);
//
//             side_fe->get_normals();
//
//             // FIXME: this is a hack to support normals at quadrature
//             // points; we don't support normals elsewhere.
//             side_fe->get_xyz();
//           }
//    }



    bool read(std::string& str);
    void add_function(std::string str);

    double operator()( const double t,
                       const double x,
                       const double y,
                       const double z,
                       const unsigned int component ) const;

    /**
     * Clears the function.
     */
    void clear () { M_expression.clear(); }

    /**
     * Returns a new copy of the function.  The new copy should be as
     * ``deep'' as necessary to allow independent destruction and
     * simultaneous evaluations of the copies in different threads.
     */
    std::unique_ptr<libMesh::FunctionBase<double> > clone () const;


    double operator() (const Point & p,
                               const double time = 0.);
    void operator() (const Point & p,
                     const double time,
                     libMesh::DenseVector<double> & output);
    double component(unsigned int i,
                     const Point & p,
                     double time = 0.0);
	void showMe(std::ostream& ofstream = std::cout );
	int size() const;
    static grammar_Type  M_grammar;
    std::vector<std::string>    M_expression;
};

} /* namespace BeatIt */

#endif /* SRC_UTIL_SPIRITFUNCTION_HPP_ */
