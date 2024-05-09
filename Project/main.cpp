#include <iostream>
#include "Eigen/Eigen"
#include "map"
#include "dfn.hpp"
#include <iomanip>

using namespace std;
using namespace FractureLibrary;

int main()
{
    FractureMesh mesh;
    string filepath = "DFN/FR3_data.txt";
    if(!ImportFR_data(filepath,mesh)){
        return 1;
    }

    for(unsigned int i=0; i<mesh.NumFractures; i++){
        cout << "Fracture id: " << i << endl;
        for(unsigned int j=0; j<mesh.MapFractures.at(i).size(); j++){
            for(unsigned int k=0; k<3; k++){
                cout << scientific << setprecision(16) <<  mesh.MapFractures.at(i)[j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    return 0;
}
