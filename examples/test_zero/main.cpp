#include <iostream>

#include "Util/IO/io.hpp"


#include <iostream>
#include <vector>
#include <fstream>
#include <string>



int main()
{
    BeatIt::printBanner(std::cout);

    std::vector<std::string> sv1;
    std::vector<std::string> sv2;
    std::vector<double> dv1;
    std::vector<double> dv2;


    std::string r1 = "va, af, fa";
    std::string r2 = "3.0, 2, 1.333";

    BeatIt::readList(r1, sv1);
    for (auto&& p : sv1) std::cout <<  p << " ";
    std::cout << std::endl;
    BeatIt::readList(r2, dv1);
    for (auto&& p : dv1) std::cout <<  p << " ";
    std::cout << std::endl;


    return 0;
}
