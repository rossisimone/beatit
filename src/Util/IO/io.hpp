#include <iostream>
#include <vector>
#include <fstream>


namespace BeatIt
{

void printBanner( std::ostream& cout);

void saveData( double time,
		       std::vector<double>& var,
			   std::ostream& output);


} // namespace BeatIt
