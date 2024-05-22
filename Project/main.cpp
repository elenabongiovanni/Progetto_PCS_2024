#include <iostream>
#include "Eigen/Eigen"
#include "map"
#include "dfn.hpp"
#include "utils.hpp"
#include <iomanip>

using namespace std;
using namespace FractureLibrary;

int main()
{
    FractureMesh mesh;
    string filepath;
    Fracture f;
    cout << "Insert filepath: ";    //scegli tra
        //DFN/FR3_data.txt
        //DFN/FR10_data.txt
        //DFN/FR50_data.txt
        //DFN/FR82_data.txt
        //DFN/FR200_data.txt
        //DFN/FR362_data.txt
    getline(cin, filepath);

    if(!ImportFR_data(filepath,mesh)){
        return 1;
    }

    for(unsigned int i=0; i<mesh.NumFractures; i++){
        //cout << "Fracture id: " << i << endl;
        for(unsigned int j=0; j<mesh.MapFractures[i].NumVertices; j++){
            for(unsigned int k=0; k<3; k++){
                //cout << scientific << setprecision(16) << mesh.MapFractures[i].vertices[j][k] << " ";
            }
            //cout << endl;
        }
        //cout << endl;
    }


    /*for(unsigned int i=0; i<mesh.NumFractures; i++){
        cout << "il baricentro della frattura " << i << " e':\n [ ";
        for(unsigned int j=0; j<3; j++){
            cout <<mesh.MapFractures.at(i).barycentre[j]<<" ";
            //findIntersections(i, mesh);
        }
        cout << " ]"<<endl;
    }*/


    for(unsigned int i=0; i<mesh.NumFractures-1; i++){
        findIntersections(i, mesh);
        cout << endl;
    }

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
