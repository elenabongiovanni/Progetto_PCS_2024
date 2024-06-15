/*#ifndef UTILS_HPP
#define UTILS_HPP

#include "dfn.hpp"

using namespace std;
using namespace Eigen;
using namespace FractureLibrary;

namespace GeometryHelpers
{
    //double findMin(const vector<double> v);

    //double findMax(const vector<double> v);

    //bool isPointInsideFracture(const Vector3d& p, const Fracture& f);

    //Vector3d calcoloPiano(Fracture& f);

    vector<Vector3d> intersezionipoligonoretta(const Vector3d& t, const Vector3d p, Fracture& f);
    vector<unsigned int> intersezionipoligonorettaLATI(const Vector3d& t, const Vector3d p, const Fracture& f, vector<unsigned int>toRemove);

    //vector<Fracture> cuttingfractures(const Fracture& f, const Trace& t, PolygonalMesh& pm);

    VectorXd PALUSolver(const MatrixXd& a, const VectorXd& b);

    bool orderLen(const Trace &a, const Trace &b);

    bool compareFirstElement(const Vector3d& a, const Vector3d& b);

    double dist(Vector3d v1, Vector3d v2);

    //double maxDist(Fracture f);

    bool onSegment(const Vector3d& p, const Vector3d& a, const Vector3d& b);

    void intersezioniSuRetta(bool& bole, vector<Vector3d>& trace, vector<Vector3d>& s1, vector<Vector3d>& s2);
    bool sameLine(const Vector3d& retta, const Fracture& f, vector<Vector3d> coordinate);


    //bool ImportFR_data(const string &filename);

    //void findIntersections(FractureMesh &mesh);

    //void printingtraces(FractureMesh& mesh, const string& file);

    //void printingfractures(FractureMesh& mesh, const string& file);

    //vector<PolygonalMesh> newpolygon(FractureMesh& mesh);

    //void printingPolygonMesh(const vector<PolygonalMesh>& pm, const string& file);

    //void defNewTrace(Trace& t, const double& d1, const double& d2, Fracture& f1, Fracture& f2, FractureMesh& fm);

}
#endif // UTILS_HPP
