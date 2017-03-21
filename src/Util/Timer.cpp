/*
 * Timer.hpp
 *
 *  Created on: 12/apr/2015
 *      Author: srossi
 */


#include "Util/Timer.hpp"

namespace BeatIt
{


Timer::Timer () : M_start(), M_end(), M_run(), M_elapsed()
{
}

void Timer::reset()
{
    M_run = false;
    M_elapsed = duration_Type::zero();
}

void Timer::restart()
{
    M_elapsed = duration_Type::zero();
    M_run = true;
    M_start = clock_Type::now();
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
    out << "Timer: elapsed time = " << elapsed().count() << "  s." << std::endl;
}

} /* namespace AthenaVMS */

