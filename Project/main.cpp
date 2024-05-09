#include <iostream>
#include "Eigen/Eigen"
#include "map"
#include "DFN.hpp"

using namespace std;
using namespace FractureLibrary;

int main()
{
    FractureMesh mesh;
    string filepath = "FR3_data.txt";
    if(!ImportFR_data(filepath,mesh)){
        return 1;
    }


    return 0;
}
