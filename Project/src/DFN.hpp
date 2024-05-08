#ifndef DFN_HPP
#define DFN_HPP

#include <iostream>
#include "Eigen/Eigen"
#include "map"

using namespace std;
using namespace Eigen;

namespace FractureLibrary
{
// struttura mesh poligonale
struct FractureMesh
{
    unsigned int NumFractures = 0; // numero fratture
    Matrix3d CoordinatesFractures = {}; // coordinate delle fratture
    vector<unsigned int> FractureId = {}; // identificatore
    unsigned int NumVertices = 0; // rete 3-dimensionale
    vector<vector<unsigned int>> Cell3DVertices = {}; // lista dei vertici
    vector<Vector2d> Cell1DVertices = {}; // descritto da due vertici: origin e end

};

bool ImportFR_data(const string &filename,
                    FractureMesh& mesh);

}
#endif // DFN_HPP
