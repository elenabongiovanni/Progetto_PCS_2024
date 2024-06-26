#ifndef DFN_H
#define DFN_H

#include <iostream>
#include "Eigen/Eigen"
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace Eigen;

namespace FractureLibrary
{


double dist(Vector3d v1, Vector3d v2);

// singola traccia condivisa da due fratture
struct Trace
{
    Vector2i fraId = {};
    vector<Vector3d> coordTrace = {};
    double len = 0.0;
    unsigned int id = 0;
    Vector3d retta = {};
    Vector3d p = {};



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
    unsigned int numTrac = 0;
    list<unsigned int> idvertici = {};
    list<unsigned int> idlati = {};
    vector<double> plane = {};
    map<unsigned int, bool> onEdge ={};
    bool isOnEdge = false;

    //calcolo distanza massima dal baricentro

    double maxDist(){
        double r = 0.0;
        for(unsigned int i=0; i< NumVertices; i++){
            double d = dist(barycentre, vertices[i]);
            if(d>r)
                r = d;
        }
        return r;
    }

    Vector3d calcoloPiano(){
        plane.clear();
        Vector3d lato1F = {};
        Vector3d lato2F = {};
        //calcolo 2 lati
        for(unsigned int j =0; j<3; j++){
            lato1F[j] = vertices[1][j] - vertices[0][j];
            lato2F[j] = vertices[3][j] - vertices[0][j];
        }

        Vector3d planeF = lato1F.cross(lato2F);
        double d = -(planeF[0] * vertices[0][0]) - (planeF[1] * vertices[0][1]) - (planeF[2] * vertices[0][2]);
        plane.reserve(4);
        for(unsigned int i=0; i<3; i++){
            plane.push_back(planeF[i]);
        }
        plane.push_back(d);
        return planeF;
    }

};


// struttura mesh poligonale
struct FractureMesh
{
    unsigned int NumFractures = 0; // numero fratture
    vector<unsigned int> FractureId = {}; // identificatore
    vector<Fracture> MapFractures = {};
    map<unsigned int, Trace> MapTrace = {};


    void printingtraces(const string& file){
        string filepath = "traces_" + file.substr(4,file.length()-4);
        ofstream outfile(filepath);
        if (outfile.is_open()) {
            outfile << "#N Traces" << endl << MapTrace.size() << endl;
            outfile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
            for(const auto &t: MapTrace){
                outfile << t.first << "; " << t.second.fraId[0] << "; " << t.second.fraId[1] << "; ";
                for(unsigned int i=0; i<2; i++){
                    for(unsigned int j=0; j<3; j++){
                        outfile << scientific << setprecision(16)<< t.second.coordTrace[i][j] << "; ";

                    }
                }
                outfile << endl;
            }
            outfile.close();
        } else {
            cout << "Impossibile aprire il file." << endl;
        }
    }

    void printingfractures(const string& file){
        string filepath = "fractures_" + file.substr(4,file.length()-4);
        ofstream outfile(filepath);
        if (outfile.is_open()) {
            for(const auto &f: MapFractures){
                outfile << "# FractureId; NumTraces" << endl << f.id << "; " << f.numTrac << endl << "# TraceId; Tips; Lenght" << endl;
                for(const auto &t: f.listPas){
                    outfile << t.id << "; false; " << scientific << setprecision(16)<< t.len << endl;
                }
                for(const auto &t: f.listNonpas){
                    outfile << t.id << "; true; " << scientific << setprecision(16)<< t.len << endl;
                }
            }
            outfile.close();
        } else {
            cout << "Impossibile aprire il file." << endl;
        }
    }


};

struct Cell0d{
    unsigned int id;
    Vector3d coordinates;
    bool old = false;
    vector<unsigned int> touched2D = {};
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
    unsigned int numFrac;

    void convertFracture(const Fracture& f, unsigned int& idVert, unsigned int& idEdge, unsigned int& id2d){
        numFrac = f.id;
        numVert = f.NumVertices;
        Cell2DVertices.resize(numVert);
        Cell2DEdges.resize(numVert);
        id = id2d++;
        for(unsigned int i=0; i<numVert; i++){
            Cell0d newCell0d(idVert,f.vertices[i]);
            newCell0d.touched2D.push_back(this->id);
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
    vector<unsigned int> Cell0DId = {};
    map<unsigned int, Cell0d> MapCell0D = {};

    unsigned int numCell1D = 0;
    vector<unsigned int> Cell1DId = {};
    map<unsigned int, Cell1d> MapCell1D = {};

    unsigned int numCell2D = 0;
    vector<unsigned int> Cell2DId = {};

    map<unsigned int, Cell2d> MapCell2D = {};
    map<unsigned int, vector<unsigned int>> MapCell2DVertices = {};
    map<unsigned int, vector<unsigned int>> MapCell2DEdges = {};
  
    void addFirstCell2d(Cell2d& c2){

        numCell0D = c2.numVert;
        numCell1D = c2.numVert;
        numCell2D = 1;
        for(Cell0d& c0: c2.Cell2DVertices){
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

        for(const Cell0d& c0d: c2.Cell2DVertices)
            MapCell2DVertices[c2.id].push_back(c0d.id);

        for(const Cell1d& c1d: c2.Cell2DEdges)
            MapCell2DEdges[c2.id].push_back(c1d.id);
        numFrac = c2.numFrac;
    }

    void addingStuff( vector<unsigned int> idLatitagliati, vector<Cell1d>& forming, vector<Cell0d>& forming0d, Cell2d& c2new1, Cell2d& c2new2, unsigned int& id2d, list<Cell2d>& next){
        unsigned int pos = 0;
        for(const unsigned int& id: idLatitagliati){
            MapCell1D.at(id).old = true;
            MapCell1D.at(id).tobecome.push_back(forming[pos++]);
            MapCell1D.at(id).tobecome.push_back(forming[pos++]);
        }

        Cell2DId.push_back(c2new1.id);
        Cell2DId.push_back(c2new2.id);
        MapCell2D[c2new1.id] = c2new1;
        MapCell2D[c2new2.id] = c2new2;

        for(const Cell0d& c0d: c2new1.Cell2DVertices)
            MapCell2DVertices[c2new1.id].push_back(c0d.id);
        for(const Cell0d& c0d: c2new2.Cell2DVertices)
            MapCell2DVertices[c2new2.id].push_back(c0d.id);
        for(const Cell1d& c1d: c2new1.Cell2DEdges)
            MapCell2DEdges[c2new1.id].push_back(c1d.id);
        for(const Cell1d& c1d: c2new2.Cell2DEdges)
            MapCell2DEdges[c2new2.id].push_back(c1d.id);

        for(Cell1d& newc1d: forming){
            Cell1DId.push_back(newc1d.id);
            MapCell1D[newc1d.id] = newc1d;
        }
        for(Cell0d& newc0d: forming0d){
            Cell0DId.push_back(newc0d.id);
            newc0d.old = true;
            MapCell0D[newc0d.id] = newc0d;
        }
        for(const unsigned int& id: idLatitagliati){
            for(unsigned int& c2: MapCell1D.at(id).touched2D){
                if(!MapCell2D.at(c2).old){

                    next.push_back(MapCell2D.at(c2));

                }
            }
        }
    }



};

bool cuttingfractures(Cell2d& f, const Trace& t, PolygonalMesh& polyMesh, list<Cell2d>& next);

bool intersLato(const Trace& t, const Cell0d& c1, const Cell0d& c2, Vector3d& inters, const PolygonalMesh& pm);

bool ImportFR_data(const string &filename, FractureMesh& mesh);

void findIntersections(FractureMesh &mesh);

void defNewTrace(Trace& t, const double& d1, const double& d2, Fracture& f1, Fracture& f2, FractureMesh& fm);

vector<Vector3d> intersezionipoligonoretta(const Vector3d& t, const Vector3d &p, vector<Vector3d> &f, bool &onEdge);

VectorXd PALUSolver(const MatrixXd& a, const VectorXd& b);

bool orderLen(const Trace &a, const Trace &b);

bool compareFirstElement(const Vector3d& a, const Vector3d& b);

bool onSegment(const Vector3d& p, const Vector3d& a, const Vector3d& b);

void intersezioniSuRetta(bool& bole, vector<Vector3d>& trace, const vector<Vector3d> &s1, const vector<Vector3d> &s2);

bool sameLine(const Vector3d& retta, const Vector3d& p, const vector<Vector3d>& f, vector<Vector3d> &coordinate);

vector<PolygonalMesh> newpolygon(FractureMesh& mesh);

void printingPolygonMesh(const vector<PolygonalMesh>& pm, const string& file);

bool checkIsNew(const vector<unsigned int>& c2d, const Vector3d& point, const PolygonalMesh& pm, unsigned int& id);

void addNewVertAndEdg(bool& firstCell2d, unsigned int& wheretoinsert, Cell2d& c2new1, Cell2d& c2new2, bool& beenFalse,
                      vector<Cell1d>& forming,vector<Cell0d>& forming0d, Vector3d& intersection, Cell1d& c1d,
                      unsigned int& vert0, unsigned int& vert1, Cell2d& f);

void addNewEdg(bool& firstCell2d, unsigned int& wheretoinsert, Cell2d& c2new1, Cell2d& c2new2, bool& beenFalse, Cell1d& c1d,
               Cell2d& f, PolygonalMesh& polyMesh, unsigned int& vert0, unsigned int& idSame);

void dividingExistingVert(const unsigned int& idSame, Cell2d& f, unsigned int& vert0, bool& firstCell2d, bool& beenFalse,
                          Cell2d& c2new1, Cell2d& c2new2, Cell1d& c1d, unsigned int& wheretoinsert );

bool cuttedByNonPas(const vector<Vector3d>& copiacoordiTrace, const Cell2d& cc, const Fracture& f, const Trace &trace);

void splitOneEdg(unsigned int& id2D, Cell2d& toCut, PolygonalMesh& pm);

}
#endif // DFN_HPP

