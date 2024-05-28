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
    unsigned int id;
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
    map<unsigned int, list<unsigned int>> MapCell1D = {};

    unsigned int numCell2D = 0;
    list<unsigned int> Cell2DId = {};
    map<unsigned int, list<unsigned int>> MapCell2DVertices = {};
    map<unsigned int, list<unsigned int>> MapCell2DEdges = {};
};


}
#endif // DFN_HPP
