#ifndef DFN_H
#define DFN_H

#include <iostream>
#include "Eigen/Eigen"
//#include "poligonalMesh.hpp"
#include <map>
#include <vector>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace Eigen;
//using namespace Polygon;

namespace FractureLibrary
{

double dist(Vector3d v1, Vector3d v2); //fatto test

// singola traccia condivisa da due fratture
struct Trace
{
    Vector2i fraId = {};
    //map<unsigned int, bool> fracturesTrace ={};
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
    unsigned int numFrac = 0;
    //map<unsigned int, vector<Vector3d>> intersectionRettaTraccia = {};
    //Vector3d normalePiano = {};
    list<unsigned int> idvertici = {};
    list<unsigned int> idlati = {};
    vector<double> plane = {};
    map<unsigned int, bool> onEdge ={};
    bool isOnEdge = false;

    //calcolo distanza massima dal baricentro
    double maxDist(){ // fatto test
        double r = 0.0;
        for(unsigned int i=0; i< NumVertices; i++){
            double d = dist(barycentre, vertices[i]);
            if(d>r)
                r = d;
        }
        return r;
    }

    Vector3d calcoloPiano(){ //fatto test
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
        cout << "d: " << d << endl;
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
                outfile << "# FractureId; NumTraces" << endl << f.id << "; " << f.numFrac << endl << "# TraceId; Tips; Lenght" << endl;
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
    unsigned int numFrac;

    void convertFracture(const Fracture& f, unsigned int& idVert, unsigned int& idEdge, unsigned int& id2d){ // fatto test
        numFrac = f.id;
        numVert = f.NumVertices;
        Cell2DVertices.resize(numVert);
        Cell2DEdges.resize(numVert);
        id = id2d++;
        //unsigned int idVert = 0;
        //unsigned int idEdge = 0;
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
    list<unsigned int> Cell0DId = {};
    map<unsigned int, Cell0d> MapCell0D = {};

    unsigned int numCell1D = 0;
    list<unsigned int> Cell1DId = {};
    map<unsigned int, Cell1d> MapCell1D = {};

    unsigned int numCell2D = 0;
    list<unsigned int> Cell2DId = {};
    map<unsigned int, Cell2d> MapCell2D = {};
    map<unsigned int, vector<unsigned int>> MapCell2DVertices = {};
    map<unsigned int, vector<unsigned int>> MapCell2DEdges = {};

    void addFirstCell2d(Cell2d& c2){ // fatto test
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

    void addingStuff(vector<bool>& angolo, vector<unsigned int> idLatitagliati, vector<Cell1d>& forming, vector<Cell0d>& forming0d, Cell2d& c2new1, Cell2d& c2new2, unsigned int& id2d, list<Cell2d>& next){
        unsigned int pos = 0;
        if(!angolo[0] && !angolo[1]){
            for(const unsigned int& id: idLatitagliati){
                MapCell1D.at(id).old = true;
                MapCell1D.at(id).tobecome.push_back(forming[pos++]);
                MapCell1D.at(id).tobecome.push_back(forming[pos++]);
            }
        }
        else if((angolo[0] && !angolo[1]) || (angolo[1] && !angolo[0])){
            MapCell1D.at(idLatitagliati[0]).old = true;
            MapCell1D.at(idLatitagliati[0]).tobecome.push_back(forming[pos++]);
            MapCell1D.at(idLatitagliati[0]).tobecome.push_back(forming[pos++]);
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
                    MapCell2D.at(c2).old = true;
                    Cell2d c2new3;
                    c2new3.id = id2d;
                    vector<unsigned int> listaidvertici;
                    vector<unsigned int> listaidlati;
                    unsigned int vert = 0;
                    for(Cell1d& lati: MapCell2D.at(c2).Cell2DEdges){
                        unsigned int idvert0 = MapCell2D.at(c2).Cell2DVertices[vert].id;
                        //unsigned int idvert1 = MapCell2D.at(c2).Cell2DVertices[(vert+1)];
                        if(MapCell1D.at(lati.id).old){
                            c2new3.Cell2DVertices.push_back(MapCell0D.at(idvert0));
                            listaidvertici.push_back(idvert0);
                            //bool toAdd=true;
                            Cell1d& c1 = MapCell1D.at(id).tobecome[0];
                            Cell1d& c2 = MapCell1D.at(id).tobecome[1];
                            if(MapCell0D.at(c1.extremes[0]).id==idvert0){
                                c2new3.Cell2DVertices.push_back(MapCell0D.at(c1.extremes[1]));
                                //c2new3.Cell2DVertices.push_back(MapCell0D.at(c2.extremes[1]));
                                c2new3.Cell2DEdges.push_back(c1);
                                listaidlati.push_back(c1.id);
                                c2new3.Cell2DEdges.push_back(c2);
                                listaidlati.push_back(c2.id);

                                /*}else if(MapCell0D.at(c1.extremes[1]).id==idvert0){
                                c2new3.Cell2DVertices.push_back(MapCell0D.at(c1.extremes[0]));
                                c2new3.Cell2DVertices.push_back(MapCell0D.at(c2.extremes[0]));
                            }else if(MapCell0D.at(c2.extremes[0]).id==idvert0){
                                c2new3.Cell2DVertices.push_back(MapCell0D.at(c2.extremes[1]));
                                c2new3.Cell2DVertices.push_back(MapCell0D.at(c1.extremes[1])); */
                            }else if(MapCell0D.at(c2.extremes[1]).id==idvert0){
                                c2new3.Cell2DVertices.push_back(MapCell0D.at(c2.extremes[0]));
                                //c2new3.Cell2DVertices.push_back(MapCell0D.at(c1.extremes[0]));
                                c2new3.Cell2DEdges.push_back(c2);
                                listaidlati.push_back(c2.id);
                                c2new3.Cell2DEdges.push_back(c1);
                                listaidlati.push_back(c1.id);
                                listaidvertici.push_back(c2.extremes[0]);
                            }

                            /*for(Cell1d& c1: MapCell1D.at(id).tobecome){
                                if(MapCell0D.at(c1.extremes[1]).id!=idvert0 && toAdd){
                                    c2new3.Cell2DVertices.push_back(MapCell0D.at(c1.extremes[1]));
                                    listaidvertici.push_back(c1.extremes[1]);
                                    toAdd=false;
                                }
                                else if(MapCell0D.at(c1.extremes[1]).id==idvert0 && toAdd){
                                    c2new3.Cell2DVertices.push_back(MapCell0D.at(c1.extremes[0]));
                                    listaidvertici.push_back(c1.extremes[0]);
                                    toAdd=false;
                                }*/

                            //c2new3.Cell2DVertices.push_back(MapCell0D.at(c1.extremes[1]));
                            /*c2new3.Cell2DEdges.push_back(c2);
                            listaidlati.push_back(c2.id);
                            c2new3.Cell2DEdges.push_back(c1);
                            listaidlati.push_back(c1.id);*/
                            //listaidvertici.push_back(c1.extremes[0]);
                            //listaidvertici.push_back(c1.extremes[1]);


                        }else{
                            c2new3.Cell2DVertices.push_back(MapCell0D.at(idvert0));
                            c2new3.Cell2DEdges.push_back(MapCell1D.at(lati.id));
                            listaidvertici.push_back(idvert0);
                            listaidlati.push_back(lati.id);
                        }
                        vert++;
                    }

                    Cell2DId.push_back(id2d);
                    c2new3.numVert = listaidvertici.size();
                    MapCell2D[id2d] = c2new3;
                    MapCell2DVertices[id2d] = listaidvertici;
                    MapCell2DEdges[id2d++] = listaidlati;
                    next.push_back(c2new3);
                }
            }
        }

    }



};

bool cuttingfractures(Cell2d& f, const Trace& t, PolygonalMesh& polyMesh, vector<Cell2d>next);

bool intersLato(const Trace& t, const Cell0d& c1, const Cell0d& c2, Vector3d& inters, const PolygonalMesh& pm); //fatto test

bool ImportFR_data(const string &filename, FractureMesh& mesh); //fatto

void findIntersections(FractureMesh &mesh);

void addingStuff(vector<bool>& angolo, vector<unsigned int> idLatitagliati, vector<Cell1d>& forming, vector<Cell0d>& forming0d, Cell2d& c2new1, Cell2d& c2new2, unsigned int& id2d, list<Cell2d>& next);

void defNewTrace(Trace& t, const double& d1, const double& d2, Fracture& f1, Fracture& f2, FractureMesh& fm); //fatto test

void printingtraces(FractureMesh& mesh, const string& file);

//void printingfractures(FractureMesh& mesh, const string& file);

vector<Vector3d> intersezionipoligonoretta(const Vector3d& t, const Vector3d &p, vector<Vector3d> &f, bool onEdge); //fatto test

vector<unsigned int> intersezionipoligonorettaLATI(const Vector3d& t, const Vector3d p, const Fracture& f, vector<unsigned int>toRemove);

//vector<Fracture> cuttingfractures(const Fracture& f, const Trace& t, PolygonalMesh& pm);

VectorXd PALUSolver(const MatrixXd& a, const VectorXd& b); //fatto

bool orderLen(const Trace &a, const Trace &b); // fatto

bool compareFirstElement(const Vector3d& a, const Vector3d& b); //fatto

bool onSegment(const Vector3d& p, const Vector3d& a, const Vector3d& b); // fatto

void intersezioniSuRetta(bool& bole, vector<Vector3d>& trace, vector<Vector3d>& s1, vector<Vector3d>& s2); //fatto test

bool sameLine(const Vector3d& retta, const Vector3d& p, const vector<Vector3d>& f, vector<Vector3d> coordinate); // fatto

vector<PolygonalMesh> newpolygon(FractureMesh& mesh);

void printingPolygonMesh(const vector<PolygonalMesh>& pm, const string& file);

bool checkIsNew(const Cell1d& c2d, const Vector3d& point, const PolygonalMesh& pm, unsigned int& id); //fatto


}
#endif // DFN_HPP
