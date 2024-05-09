#ifndef DFN_H
#define DFN_H

#include <iostream>
#include "Eigen/Eigen"
#include <map>
#include <vector>

using namespace std;
using namespace Eigen;

namespace FractureLibrary
{
// struttura mesh poligonale
struct FractureMesh
{
    unsigned int NumFractures = 0; // numero fratture
    //Matrix3d CoordinatesFractures = {}; // coordinate delle fratture
    vector<unsigned int> FractureId = {}; // identificatore
    unsigned int NumVertices = 0; // rete 3-dimensionale
    //vector<vector<unsigned int>> Cell3DVertices = {}; // lista dei vertici
    //vector<Vector2d> Cell1DVertices = {}; // descritto da due vertici: origin e end
    map<unsigned int, vector<Vector3d>> MapFractures = {};
    //vector<Vector3d<double>> CoordinatesFractures = {};

};

bool ImportFR_data(const string &filename,
                   FractureMesh& mesh);

}
#endif // DFN_HPP
