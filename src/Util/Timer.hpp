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
    typedef std::chrono::milliseconds milliseconds;
    typedef std::chrono::duration<double> duration_Type;
public:

    Timer ();

    void reset();
    void start();
    void stop();
    duration_Type elapsed();
    void print( std::ostream& out);

    duration_Type  M_elapsed;
    timePoint_Type M_start;
    timePoint_Type M_end;
    bool           M_run;
};



Timer::Timer () : M_start(), M_end(), M_run(), M_elapsed()
{
}

void Timer::reset()
{
    M_run = false;
    M_elapsed = duration_Type::zero();
}

void Timer::start()
{
    M_run = true;
    M_start = clock_Type::now();
}

void Timer::stop()
{
	if(M_run)
	{
		M_run = false;
	    M_end = clock_Type::now();
	    M_elapsed += M_end - M_start;
	}
}


Timer::duration_Type
Timer::elapsed()
{
    if( !M_run) return M_elapsed;
    else
    {
        duration_Type elapsed = M_elapsed;
        auto end = clock_Type::now();
        elapsed += end - M_start;
        return elapsed;
    }
}

void Timer::print( std::ostream& out)
{
    out << "Timer: elapsed time = " << elapsed().count() << " ms." << std::endl;
}

} /* namespace AthenaVMS */

#endif /* TESTSUITE_TEST_SMALLOPERATIONS_TEMPORARIES_TIMER_HPP_ */
