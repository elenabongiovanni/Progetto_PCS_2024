#ifndef UTILS_HPP
#define UTILS_HPP

#include "dfn.hpp"

using namespace std;
using namespace Eigen;

namespace FractureLibrary
{

bool ImportFR_data(const string &filename,
                   FractureMesh& mesh);

double dist(Vector3d v1, Vector3d v2);

double maxDist(Fracture f);

void findIntersections(const unsigned int &id, FractureMesh &mesh);



}
#endif // UTILS_HPP
