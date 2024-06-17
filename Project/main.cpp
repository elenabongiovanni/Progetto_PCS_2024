#include <iostream>
#include "Eigen/Eigen"
#include "map"
#include "DFN.hpp"
//#include "utils.hpp"
#include <iomanip>
#include <fstream>
//#include "poligonalMesh.hpp"

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
    string filepath = "FR10_data.txt";


    if(!ImportFR_data(filepath, mesh)){
        return 2;
    }

    /*for(unsigned int i=0; i<mesh.NumFractures; i++){
        cout << "Fracture id: " << i << endl;
        for(unsigned int j=0; j<mesh.MapFractures[i].NumVertices; j++){
            for(unsigned int k=0; k<3; k++){
                cout << scientific << setprecision(16) << mesh.MapFractures[i].vertices[j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }*/


    /*for(unsigned int i=0; i<mesh.NumFractures; i++){
        cout << "il baricentro della frattura " << i << " e':\n [ ";
        for(unsigned int j=0; j<3; j++){
            cout <<mesh.MapFractures.at(i).barycentre[j]<<" ";
            //findIntersections(i, mesh);
        }
        cout << " ]"<<endl;
    }*/


    //for(unsigned int i=0; i<mesh.NumFractures-1; i++){
    findIntersections(mesh);
    mesh.printingtraces(filepath);
    mesh.printingfractures(filepath);

    vector<PolygonalMesh> pp = newpolygon(mesh);
    printingPolygonMesh(pp, filepath);

    /*for(unsigned int id: p.Cell0DId){
        cout << id << " ";
    }
    cout << endl;
    cout<< p.numCell0D << endl;
    cout<< p.numCell1D << endl;
    cout<< p.numCell2D << endl;
        //cout << endl;
    //}

    /*for(unsigned int i=0; i<mesh.vecTrace.size(); i++){
        Trace t = mesh.vecTrace[i];

        for(auto pair : t.fracturesTrace){
            cout << "bau" << endl;
            if(pair.second){
                cout <<"fratturna numero "<< pair.first << " e' passante" << endl;
            }
        }
    }*/


    return 0;
}
