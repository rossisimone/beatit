/*
 * Timer.hpp
 *
 *  Created on: 12/apr/2015
 *      Author: srossi
 */

#ifndef TIMER_HPP_
#define TIMER_HPP_


#include <chrono>
#include <ostream>


namespace BeatIt
{

class Timer
{
    typedef std::chrono::high_resolution_clock::time_point timePoint_Type;
    typedef std::chrono::high_resolution_clock clock_Type;
    typedef std::chrono::seconds seconds;
    typedef std::chrono::duration<double > duration_Type;
public:

    Timer ();

    void reset();
    void restart();
    void start();
    void stop();
    duration_Type elapsed();
    void print( std::ostream& out);

    duration_Type  M_elapsed;
    timePoint_Type M_start;
    timePoint_Type M_end;
    bool           M_run;
};



} /* namespace AthenaVMS */

#endif /* TESTSUITE_TEST_SMALLOPERATIONS_TEMPORARIES_TIMER_HPP_ */
