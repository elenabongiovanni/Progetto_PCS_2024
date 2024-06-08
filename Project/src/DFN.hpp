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

// singola traccia condivisa da due fratture
struct Trace
{
    Vector2i fraId = {};
    //map<unsigned int, bool> fracturesTrace ={};
    vector<Vector3d> coordTrace = {};
    double len = 0.0;
    unsigned int id = 0;
    Vector3d retta = {};

};


// singola frattura
struct Fracture
{
    unsigned int id = 0;
    unsigned int NumVertices = 0;
    vector<Vector3d> vertices = {};
    Vector3d barycentre = {};
    list<Trace> listPas = {};
    list<Trace> listNonpas = {};
    map<unsigned int, bool> tips = {};
    unsigned int numFrac = 0;
    //map<unsigned int, vector<Vector3d>> intersectionRettaTraccia = {};
    //Vector3d normalePiano = {};
    list<unsigned int> idvertici = {};
    list<unsigned int> idlati = {};
    vector<double> plane = {};
    map<unsigned int, bool> onEdge ={};
    bool isOnEdge = false;

};


// struttura mesh poligonale
struct FractureMesh
{
    unsigned int NumFractures = 0; // numero fratture
    //Matrix3d CoordinatesFractures = {}; // coordinate delle fratture
    vector<unsigned int> FractureId = {}; // identificatore
    //unsigned int NumVertices = 0; // rete 3-dimensionale //il numero di vertici dovrebe essere un elemento della singola frattura non della mesh
    //vector<vector<unsigned int>> Cell3DVertices = {}; // lista dei vertici
    //vector<Vector2d> Cell1DVertices = {}; // descritto da due vertici: origin e end
    vector<Fracture> MapFractures = {}; // potresti cambiare con freactures
    //vector<Fracture> MapFractures = {};
    //vector<Vector3d<double>> CoordinatesFractures = {};
    //vector<Trace> vecTrace = {};
    map<unsigned int, Trace> MapTrace = {};


};

struct PolygonalMesh
{
    unsigned int numFrac;

    unsigned int numCell0D = 0;
    list<unsigned int> Cell0DId = {};
    map<unsigned int, Vector3d> MapCell0D = {};

    unsigned int numCell1D = 0;
    list<unsigned int> Cell1DId = {};
    map<unsigned int, vector<unsigned int>> MapCell1D = {};

    unsigned int numCell2D = 0;
    list<unsigned int> Cell2DId = {};
    map<unsigned int, list<unsigned int>> MapCell2DVertices = {};
    map<unsigned int, list<unsigned int>> MapCell2DEdges = {};
};


vector<Fracture> cuttingfractures(Fracture& f, const Trace& t, PolygonalMesh& pm);

bool ImportFR_data(const string &filename, FractureMesh& mesh);

void findIntersections(FractureMesh &mesh);

vector<PolygonalMesh> newpolygon(FractureMesh& mesh);

void printingPolygonMesh(const vector<PolygonalMesh>& pm, const string& file);

void defNewTrace(Trace& t, const double& d1, const double& d2, Fracture& f1, Fracture& f2, FractureMesh& fm);

void printingtraces(FractureMesh& mesh, const string& file);

void printingfractures(FractureMesh& mesh, const string& file);



}
#endif // DFN_HPP
