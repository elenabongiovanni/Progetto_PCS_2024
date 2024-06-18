#include <iostream>
#include "Eigen/Eigen"
#include "map"
#include "dfn.hpp"
#include <iomanip>
#include <fstream>


using namespace std;
using namespace Eigen;
using namespace FractureLibrary;
//using namespace Polygon;

int main(int argc, char ** argv)
{
    FractureMesh mesh;

    if(argc == 1)
        return 1;

    istringstream str(argv[1]);
    string name;
    str >> name;

    string filepath = "DFN/" + name + ".txt";



    if(!ImportFR_data(filepath, mesh)){
        return 2;
    }
     
    findIntersections(mesh);
    mesh.printingtraces(filepath);
    mesh.printingfractures(filepath);

    vector<PolygonalMesh> pp = newpolygon(mesh);
    printingPolygonMesh(pp, filepath);

    


    return 0;

}

