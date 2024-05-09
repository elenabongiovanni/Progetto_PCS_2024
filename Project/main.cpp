#include <iostream>
#include "Eigen/Eigen"
#include "map"
#include "DFN.hpp"

using namespace std;
using namespace FractureLibrary;

int main()
{
    FractureMesh mesh;
    string filepath = "FR3_data";
    if(!ImportFR_data(filepath,mesh)){
        return 5;
    }


    return 0;
}
