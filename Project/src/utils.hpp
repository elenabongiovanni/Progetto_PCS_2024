#ifndef UTILS_HPP
#define UTILS_HPP

#include "dfn.hpp"

using namespace std;
using namespace Eigen;

namespace FractureLibrary
{

VectorXd PALUSolver(const MatrixXd& a, const VectorXd& b);

bool compareFirstElement(const Vector3d& a, const Vector3d& b);

double dist(Vector3d v1, Vector3d v2);

double maxDist(Fracture f);

bool onSegment(const Vector3d& p, const Vector3d& a, const Vector3d& b);

bool ImportFR_data(const string &filename, FractureMesh& mesh);

void findIntersections(const unsigned int &id, FractureMesh &mesh);


}
#endif // UTILS_HPP
