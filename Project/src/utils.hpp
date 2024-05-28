#ifndef UTILS_HPP
#define UTILS_HPP

#include "dfn.hpp"

using namespace std;
using namespace Eigen;

namespace FractureLibrary
{
    //double findMin(const vector<double> v);

    //double findMax(const vector<double> v);

    bool isPointInsideFracture(const Vector3d& p, const Fracture& f);

    Vector3d calcoloPiano(const Fracture& f);

    vector<Vector3d> intersezionipoligonoretta(const Vector3d& t, const Vector3d p, const Fracture& f);

    vector<Fracture> cuttingfractures(const Fracture& f, const Trace& t);

    VectorXd PALUSolver(const MatrixXd& a, const VectorXd& b);

    bool orderLen(const Trace &a, const Trace &b);

    bool compareFirstElement(const Vector3d& a, const Vector3d& b);

    double dist(Vector3d v1, Vector3d v2);

    double maxDist(Fracture f);

    bool onSegment(const Vector3d& p, const Vector3d& a, const Vector3d& b);

    bool ImportFR_data(const string &filename, FractureMesh& mesh);

    void findIntersections(FractureMesh &mesh);

    void printingtraces(FractureMesh& mesh, const string& file);

    void printingfractures(FractureMesh& mesh, const string& file);

    PolygonalMesh newpolygon(FractureMesh& mesh);

}
#endif // UTILS_HPP
