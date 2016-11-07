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


#ifndef SRC_UTIL_GRAMMAR_HPP_
#define SRC_UTIL_GRAMMAR_HPP_


#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include "boost/function.hpp"
#include "boost/bind.hpp"
#include <boost/spirit/include/qi_symbols.hpp>
#include <string>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/phoenix/stl/cmath.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <cmath>
#include <boost/spirit/include/qi_grammar.hpp>


#define _USE_MATH_DEFINES

namespace qi = boost::spirit::qi;
namespace ascii=boost::spirit::ascii;
namespace phoenix = boost::phoenix;
using boost::spirit::ascii::space;
using boost::spirit::qi::symbols;


namespace BeatIt
{
	constexpr static double S_PI = M_PI;
	constexpr static double S_E = M_E;
}


template< typename Iterator >
struct Grammar : public virtual qi::grammar<  Iterator, double(), ascii::space_type >
{
    typedef qi::symbols<char, double >  variable_Type;


    Grammar() : Grammar::base_type(M_expression)
    {
        M_symbol.add("pi", BeatIt::S_PI);
        M_symbol.add("e",  BeatIt::S_E);

        using qi::double_;
        using qi::_val;
        using qi::_1;


//        start =
//                ( assign
//                    | M_expression );
//        M_expression =
//                *compare                                    [qi::_val = qi::_1]
//             |
//                *plus_minus        [qi::_val = qi::_1]                          [qi::_val = qi::_1]
//            ;


//        M_expression =
//                *compare                                    [qi::_val = qi::_1]
//                |
//                (
//                 plus_minus                       [_val = _1]
//                )
//            >> *(   ('+' >> plus_minus            [_val += _1])
//                |   ('-' >> plus_minus            [_val -= _1])
//                )
//                |
//                *compare                                    [qi::_val = qi::_1]
//            ;
//
//
//        compare =
//                plus_minus                              [qi::_val = qi::_1]
//            >> * (
//                    qi::lit (">=") >> plus_minus     [qi::_val = qi::_val >= qi::_1]
//                |   qi::lit ("<=") >> plus_minus     [qi::_val = qi::_val <= qi::_1]
//                |   qi::lit (">")  >> plus_minus     [qi::_val = qi::_val >  qi::_1]
//                |   qi::lit ("<")  >> plus_minus     [qi::_val = qi::_val <  qi::_1]
//            )
//            ;
//
//
//        plus_minus =
//                (
//                    multiply_divide                  [_val = _1]
//                )
//            >> *(   ('*' >> multiply_divide          [_val *= _1])
//                |   ('/' >> multiply_divide          [_val /= _1])
//                )
//            ;
//
//        multiply_divide =
//                 (
//                           '-' >> power   [_val = _1]
//                    >>  (
//                           '^' >> power   [_val = -phoenix::pow(_val, _1)]
//                        )
//                    >> * (
//                           '^' >> power   [_val = phoenix::pow(_val, _1)] //[qi::_val = -phoenix::bind (&Grammar::pow, this, qi::_val, qi::_1)]
//
//                     )
//                 )
//                 |
//                 (
//                            power         [_val = _1]
//            >> * (
//                    '^' >>  power         [_val = phoenix::pow(_val, _1)]
//                 )
//
//                 )
//            ;
//
//
//        power =
//            qi::double_                     [_val = _1]
//            | function                      [_val = _1]
//            | group                         [_val = _1]
//            |   ('-' >> power               [_val =-_1])
//            |   ('+' >> power               [_val = _1])
//            | M_symbol                      [_val = _1]
//            ;
//


        M_expression =
            *compare                                [qi::_val = qi::_1]
            ;


        compare =
                plus_minus                              [qi::_val = qi::_1]
            >> * (
                    qi::lit (">=") >> plus_minus     [qi::_val = qi::_val >= qi::_1]
                |   qi::lit ("<=") >> plus_minus     [qi::_val = qi::_val <= qi::_1]
                |   qi::lit (">")  >> plus_minus      [qi::_val = qi::_val > qi::_1]
                |   qi::lit ("<")  >> plus_minus      [qi::_val = qi::_val < qi::_1]
                 );

        plus_minus =
                multiply_divide                         [qi::_val = qi::_1]
            >> * (
                qi::lit ('+') >> multiply_divide [qi::_val += qi::_1]
                |
                qi::lit ('-') >> multiply_divide [qi::_val -= qi::_1]
            )
            ;


        multiply_divide =
                power                                [qi::_val = qi::_1]
            >> * (
                qi::lit ('*') >> power        [qi::_val *= qi::_1]
                |
                qi::lit ('/') >> power        [qi::_val /= qi::_1]
            )
            ;


        power =
            (
                part                            [qi::_val = qi::_1]
                >> * (
                    qi::lit ('^') >> part      [_val = phoenix::pow(_val, _1)]
                )
            )
            ;
        part =
                    qi::double_                     [_val = _1]
                    | function                      [_val = _1]
                    | group                         [_val = _1]
                    |   ('-' >> power               [_val =-_1])
                    |   ('+' >> power               [_val = _1])
                    | M_symbol                      [_val = _1]
                    ;
        function = (
                       "abs"  >> group      [_val = phoenix::fabs(_1)]
                   |
                       "cos"  >> group      [_val = phoenix::cos(_1)]
                   |
                       "sin"  >> group      [_val = phoenix::sin(_1)]
                   |
                       "tan"  >> group      [_val = phoenix::tan(_1)]
                   |
                       "cosh" >> group      [_val = phoenix::cosh(_1)]
                   |
                       "sinh" >> group      [_val = phoenix::sinh(_1)]
                   |
                       "tanh" >> group      [_val = phoenix::tanh(_1)]
                   |
                       "acos" >> group      [_val = phoenix::acos(_1)]
                   |
                       "asin" >> group      [_val = phoenix::asin(_1)]
                   |
                       "atan" >> group      [_val = phoenix::atan(_1)]
                   |
                       "sqrt" >> group      [_val = phoenix::sqrt(_1)]
                   |
                       "log"  >> group      [_val = phoenix::log(_1)]
                   |
                       "log10">> group      [_val = phoenix::log10(_1)]
                   |
                       "exp"  >> group      [_val = phoenix::exp(_1)]
                   );

        group = '(' >> M_expression   [_val = _1] >> ')' ;
  }
    

    void setSymbol(std::string name, double value)
    {
        double* p = M_symbol.find ( name );
        if ( p != 0 )
        {
            *p = value;
        }
        else
        {
            M_symbol.add ( name, value );
        }
    }



    qi::rule<Iterator, double(), ascii::space_type> M_expression;
    qi::rule<Iterator, double(), ascii::space_type> multiply_divide;
    qi::rule<Iterator, double(), ascii::space_type> plus_minus;
    qi::rule<Iterator, double(), ascii::space_type> power;
    qi::rule<Iterator, double(), ascii::space_type> part;
    qi::rule<Iterator, double(), ascii::space_type> function;
    qi::rule<Iterator, double(), ascii::space_type> group;
    qi::rule<Iterator, double(), ascii::space_type> compare;
    qi::rule<Iterator, void(),   ascii::space_type > assign;
    qi::rule<Iterator, double(), ascii::space_type> start;

    variable_Type M_symbol;
};

#endif /* SRC_UTIL_GRAMMAR_HPP_ */

