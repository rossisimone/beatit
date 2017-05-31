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
 * \file SpiritFunction.cpp
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

#include "Util/SpiritFunction.hpp"
#include "libmesh/point.h"
#include "Util/IO/io.hpp"

namespace BeatIt
{

SpiritFunction::grammar_Type   SpiritFunction::M_grammar;

SpiritFunction::SpiritFunction()
    : M_expression()
{
    _initialized = true;
    _is_time_dependent = true;
}

SpiritFunction::~SpiritFunction()
{

}

bool
SpiritFunction::read(std::string& str)
{
    M_expression.clear();

//    auto f( str.begin() );
//    auto l( str.end()   );
//    namespace qi = boost::spirit::qi;
//    namespace ascii = boost::spirit::ascii;
//    bool ok = qi::phrase_parse( f, l,  ( *~boost::spirit::qi::char_(",") ) % ',', ascii::blank, M_expression);
    bool ok = BeatIt::readList(str, M_expression);
    _initialized = ok;
    if (!ok )
        throw std::runtime_error("parser error: '" + str ); // FIXME
    return ok;
}

void
SpiritFunction::showMe(std::ostream& ofstream)
{
	ofstream << "\t\t\t size: " << M_expression.size() << std::endl;
	for(auto&& str : M_expression)	ofstream << " \t\t\t function " << str << std::endl;
}


void
SpiritFunction::add_function(std::string str)
{
    _initialized = true;
    M_expression.push_back(str);
}


double
SpiritFunction::operator()(  const double t,
                            const double x,
                            const double y,
                            const double z,
                            const unsigned int component ) const
{

    double result = 0.0;

    if( component < M_expression.size() )
    {
        M_grammar.setSymbol("t", t);
        M_grammar.setSymbol("x", x);
        M_grammar.setSymbol("y", y);
        M_grammar.setSymbol("z", z);

        iterator iter = M_expression[component].begin();
        iterator end = M_expression[component].end();
        bool r = phrase_parse(iter, end, M_grammar, space, result);
        if(!r)
        {
            std::cout << "\nSpiritFunction::Parse failed!!!!" << std::endl;
            exit(-1);
        }
    }
    else
    {
        std::cout << "SpiritFunction: asked for component " << component
        		         << " but I have only " <<   M_expression.size() << " components." << std::endl;
    }
    return result;

}


libMesh::UniquePtr<libMesh::FunctionBase<double> >
SpiritFunction::clone () const
{
//    libMesh::UniquePtr< SpiritFunction > pcopy;
//    pcopy.reset(new SpiritFunction() );
//    pcopy->_master = _master;
//    pcopy->_is_time_dependent = true;
//    pcopy->_initialized = true;
//    for(int i = 0; i < M_expression.size(); ++i )
//    pcopy->M_expression.push_back(M_expression[i]);
//    return static_cast<libMesh::UniquePtr<libMesh::FunctionBase<double> > >(pcopy);
    SpiritFunction* pcopy = new SpiritFunction();
    pcopy->_master = _master;
    pcopy->_is_time_dependent = true;
    pcopy->_initialized = true;
    for(int i = 0; i < M_expression.size(); ++i )
    {
        pcopy->M_expression.push_back(M_expression[i]);
    }
    return libMesh::UniquePtr<libMesh::FunctionBase<double> >(pcopy);
}


double
SpiritFunction::operator() (const Point & p,
                           const double time )
{
    return component(0, p, time );
}



void
SpiritFunction::operator() (const Point & p,
                 const double time,
                 libMesh::DenseVector<double> & output)
{
    auto n = output.size();
    if(n != M_expression.size())
    {
        throw std::runtime_error("SpiritFunction::operator() output and M_expression have different sizes.");
    }
    int i = 0;
    for(; i <  n; i++ ) output(i) = component(i, p, time );
}

double
SpiritFunction::component(unsigned int i,
                 const Point & p,
                 double time)
{
 return operator ()(time, p(0), p(1), p(2), i);
}

int
SpiritFunction::size() const
{
	return M_expression.size();
}


} /* namespace BeatIt */
