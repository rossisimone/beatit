#include "Util/IO/io.hpp"

namespace BeatIt
{



void printBanner( std::ostream& cout)
{
	cout << "\n\t .______    _______     ___   .___________.    __  .___________.";
	cout << "\n\t |   _  \\  |   ____|   /   \\  |           |   |  | |           |";
	cout << "\n\t |  |_)  | |  |__     /  ^  \\ `---|  |----`   |  | `---|  |---- ";
	cout << "\n\t |   _  <  |   __|   /  /_\\  \\    |  |        |  |     |  |     ";
	cout << "\n\t |  |_)  | |  |____ /  _____  \\   |  |        |  |     |  |     ";
	cout << "\n\t |______/  |_______/__/     \\__\\  |__|        |__|     |__| 	 ";
	cout << "\n";
}


void saveData( double time,
		       std::vector<double>& var,
			   std::ostream& output)
{
	output << time;
	for(const auto& v : var)
	{
		output << " " << v;
	}
	output << "\n";
}


} // namespace BeatIt


