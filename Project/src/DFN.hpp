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
    //map<unsigned int, bool> fracturesTrace ={};
    vector<Vector3d> coordTrace = {};
    double length = 0.0;
    Vector2i traceId= {};

};


// singola frattura
struct Fracture
{
    unsigned int id = 0;
    unsigned int NumVertices = 0;
    vector<Vector3d> vertices = {};
    Vector3d barycentre = {};
    //map<unsigned int, vector<Vector3d>> Intersections = {};
    //metti strutture dati ordinabile a poco costo  (quicksort?)
    //traccepassanti
    //tracceNONpassanti
    list<Trace> listPas = {};
    list<Trace> listNonpas = {};
    vector <bool> Tips;

    //vector<Trace> vecPas = {};
    //vector<Trace> vecNonpas = {};

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
    map<unsigned int, Fracture> MapFractures = {}; // potresti cambiare con freactures
    //vector<Fracture> MapFractures = {};
    //vector<Vector3d<double>> CoordinatesFractures = {};
    vector<Trace> vecTrace = {};
    map<unsigned int, Trace> MapTrace = {};


};


}
#endif // DFN_HPP
