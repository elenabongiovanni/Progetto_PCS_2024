/*#ifndef POLIGONALMESH_HPP
#define POLIGONALMESH_HPP
#include <iostream>
#include "Eigen/Eigen"
#include <map>
#include <vector>
#include "dfn.hpp"
using namespace std;
using namespace Eigen;
using namespace FractureLibrary;

namespace Polygon{

struct Cell0d{
    unsigned int id;
    Vector3d coordinates;
    bool old = false;

    //Cell0d(): id(0), coordinates{} {};
    Cell0d():id(0), coordinates{} {}
    Cell0d(unsigned int& i, const Vector3d vert): id(i), coordinates(vert) {}
};

struct Cell1d{
    unsigned int id;
    vector<unsigned int> extremes = {};
    vector<unsigned int> touched2D = {};
    bool old = false;
    vector<Cell1d> tobecome = {};

    Cell1d(): id(0), extremes{}, touched2D{} {}
    Cell1d(unsigned int& i, const vector<unsigned int> vert): id(i), extremes(vert), old(false) {}
};

struct Cell2d{
    unsigned int id;
    vector<Cell0d> Cell2DVertices;
    vector<Cell1d> Cell2DEdges;
    unsigned int numVert;
    bool old = false;

    void convertFracture(const Fracture& f, unsigned int& idVert, unsigned int& idEdge){
        numVert = f.NumVertices;
        Cell2DVertices.resize(numVert);
        Cell2DEdges.resize(numVert);
        /*unsigned int idVert = 0;
        unsigned int idEdge = 0;*/
        /*for(unsigned int i=0; i<numVert; i++){
            Cell0d newCell0d(idVert,f.vertices[i]);
            idVert++;
            Cell2DVertices[i] = newCell0d;
        }
        for(unsigned int i=0; i<numVert; i++){
            unsigned int j=(i+1)%numVert;

            Cell1d newCell1d(idEdge, {Cell2DVertices[i].id, Cell2DVertices[j].id});
            idEdge++;
            newCell1d.touched2D.push_back(this->id);
            Cell2DEdges[i] = newCell1d;
        }
    }
};

struct PolygonalMesh
{
    unsigned int numFrac;

    unsigned int numCell0D = 0;
    list<unsigned int> Cell0DId = {};
    map<unsigned int, Cell0d> MapCell0D = {};

    unsigned int numCell1D = 0;
    list<unsigned int> Cell1DId = {};
    map<unsigned int, Cell1d> MapCell1D = {};

    unsigned int numCell2D = 0;
    list<unsigned int> Cell2DId = {};
    map<unsigned int, Cell2d> MapCell2D = {};
    map<unsigned int, list<unsigned int>> MapCell2DVertices = {};
    map<unsigned int, list<unsigned int>> MapCell2DEdges = {};

    void addFirstCell2d(const Cell2d& c2){
        numCell0D = c2.numVert;
        numCell1D = c2.numVert;
        numCell2D = 1;
        for(Cell0d c0: c2.Cell2DVertices){
            Cell0DId.push_back(c0.id);
            MapCell0D[c0.id] = c0;
            c0.old = true;
        }
        for(const Cell1d& c1: c2.Cell2DEdges){
            Cell1DId.push_back(c1.id);
            MapCell1D[c1.id] = c1;
        }
        Cell2DId.push_back(c2.id);
        MapCell2D[c2.id] = c2;
    }



};

    bool cuttingfractures(Cell2d& f, const Trace& t, PolygonalMesh& polyMesh, list<Cell2d> newCell, vector<Cell2d>next);
    void printingPolygonMesh(const vector<PolygonalMesh>& pm, const string& file);
    bool intersLato(const Vector3d& t, const Cell1d c1d, bool angolo, Vector3d& inters, PolygonalMesh& pm);


}
#endif // POLIGONALMESH_HPP*/
