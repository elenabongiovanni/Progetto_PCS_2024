#include "dfn.hpp"
//#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
#include <algorithm>
#include <map>
#include <iomanip>
#include <list>
#include "poligonalMesh.hpp"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
//using namespace GeometryHelpers;
//using namespace Polygon;

namespace FractureLibrary {

double tol = 10 * numeric_limits<double>::epsilon();
unsigned int idTrace = 0;

void defNewTrace(Trace& t, const double& d1, const double& d2, Fracture& f1, Fracture& f2, FractureMesh& fm){
    double lunghezzatraccia = dist(t.coordTrace[0],t.coordTrace[1]);
    t.len = lunghezzatraccia;
    t.id = idTrace;
    fm.MapTrace[idTrace++] = t;
    f1.numFrac++;
    f2.numFrac++;
    //t.retta = t;

    if(abs(lunghezzatraccia-d1)<tol){
        f1.listPas.push_back(t);
        f1.tips[t.id] = false;
    }else{
        f1.listNonpas.push_back(t);
        f1.tips[t.id] = true;
    }
    if(abs(lunghezzatraccia - d2)<tol){
        f2.listPas.push_back(t);
        f2.tips[t.id] = false;
    }else{
        f2.listNonpas.push_back(t);
        f2.tips[t.id] = true;
    }
}

bool  sameLine(const Vector3d& ret, const Vector3d& p, const vector<Vector3d>& f, vector<Vector3d> coordinate){
    double norma1 = ret.norm();
    Vector3d rettaN = {ret[0]/norma1, ret[1]/norma1, ret[2]/norma1};
    for(int k=0; k<f.size(); k++){
        int vert0 = k;
        int vert1 = (k+1)%f.size();
        Vector3d retta2 = {f[vert1][0]-f[vert0][0],
                           f[vert1][1]-f[vert0][1],
                           f[vert1][2]-f[vert0][2]};
        double norma2 = retta2.norm();
        Vector3d retta2N = {retta2[0]/norma2, retta2[1]/norma2, retta2[2]/norma2};
        if(abs(abs(rettaN[0])-abs(retta2N[0]))<tol &&
            abs(abs(rettaN[1])-abs(retta2N[1]))<tol &&
            abs(abs(rettaN[2])-abs(retta2N[2]))<tol &&
            onSegment(p, f[vert0], f[vert1])){
            coordinate.push_back(f[vert0]);
            coordinate.push_back(f[vert1]);
            return true;
        }
    }
    return false;
}




bool intersLato(const Trace& tra, const Cell0d& c1, const Cell0d& c2, bool angolo, Vector3d& inters, const PolygonalMesh& pm){
    Vector3d t = tra.retta;
    //t[0] = tra.coordTrace[1][0] - tra.coordTrace[1][0];
    //t[1] = tra.coordTrace[1][1] - tra.coordTrace[1][1];
    //t[2] = tra.coordTrace[1][2] - tra.coordTrace[1][2];
    //double norma = t.norm();
    //t[0] = t[0]/norma;
    //t[1] = t[1]/norma;
    //t[2] = t[2]/norma;
    Vector3d p = tra.coordTrace[0];
    Vector3d tr;
    Vector3d p1 = c1.coordinates;
    Vector3d p2 = c2.coordinates;
    tr[0] = p2[0] - p1[0];
    tr[1] = p2[1] - p1[1];
    tr[2] = p2[2] - p1[2];

    Vector3d prod_t_t1 = t.cross(tr);
    Vector3d prod_t_p = (p - p2).cross(t);
    //Vector3d prod_t1_p = (f.vertices[vert1]-p).cross(tr);

    double alfa = (prod_t_p.dot(prod_t_t1))/(prod_t_t1.dot(prod_t_t1));
    //double beta = (prod_t1_p.dot(prod_t_t1))/(prod_t_t1.dot(prod_t_t1));

    if(alfa>=0-tol && alfa <=1+tol){
        for(unsigned int i=0; i<3; i++){
            inters[i] = p1[i]*alfa + (1-alfa)*p2[i];
        }
        if(abs(alfa)<tol || abs(alfa-1)<tol){
            angolo = true;
        }
        return true;
    }
    else
        return false;
}

/*vector<unsigned int> intersezionipoligonorettaLATI(const Vector3d& t, const Vector3d p, const Fracture& f, vector<unsigned int>toRemove){
    vector<unsigned int> intersectionsF = {};
    intersectionsF.reserve(2);
    unsigned int index = 0;
    for(unsigned int l = 0; l < f.NumVertices; l++){
        unsigned int vert0 = l;
        unsigned int vert1 = (l+1)%f.NumVertices;

        Vector3d tr;
        tr[0] = f.vertices[vert1][0] - f.vertices[vert0][0];
        tr[1] = f.vertices[vert1][1] - f.vertices[vert0][1];
        tr[2] = f.vertices[vert1][2] - f.vertices[vert0][2];

        Vector3d prod_t_t1 = t.cross(tr);
        Vector3d prod_t_p = (p-f.vertices[vert1]).cross(t);
        //Vector3d prod_t1_p = (f.vertices[vert1]-p).cross(tr);

        double alfa = (prod_t_p.dot(prod_t_t1))/(prod_t_t1.dot(prod_t_t1));
        //double beta = (prod_t1_p.dot(prod_t_t1))/(prod_t_t1.dot(prod_t_t1));

        if(alfa>=0 && alfa <=1){

            intersectionsF.push_back(l + 1 + index);
            index++;
        }

        if(intersectionsF.size()==2){
            if(abs(intersectionsF[0][0] - intersectionsF[1][0])<tol &&  abs(intersectionsF[0][1] - intersectionsF[1][1])<tol &&
                abs(intersectionsF[0][2] - intersectionsF[1][2])<tol ){
                intersectionsF.pop_back();
                continue;
            }
            break;  //esco dal ciclo perchè di nuovo non è posisbile che una retta intersechi un poligono in più di due punti
        }           // potremmo anche mettere questo if fuori dall'else e levare il break prima , come ci piace di più
    }
    return intersectionsF;
}*/

vector<Vector3d> intersezionipoligonoretta(const Vector3d& t, const Vector3d& p, vector<Vector3d>& f, bool onEdge){
    vector<Vector3d> intersectionsF = {};
    intersectionsF.reserve(2);


    /*Vector3d plane = f.calcoloPiano();
    if(abs(plane.dot(t))<10e-3 && abs(f.plane[0]*p[0]+f.plane[1]*p[1]+f.plane[2]*p[2]+f.plane[3])<10e-3){
        cout <<"stanno sullo stesso piano."<<endl;
    }
    else{
        cout << "c'è qualquadra che non cosa" << endl;
        cout << abs(plane.dot(t)) << endl;
        cout << abs(f.plane[0]*p[0]+f.plane[1]*p[1]+f.plane[2]*p[2]+f.plane[3]) << endl;
    }*/
    //controllo che non la traccia non sia su un lato
    if(sameLine(t, p ,f,intersectionsF)){
        onEdge = true;
        return intersectionsF;
    }
    for(unsigned int l = 0; l < f.size(); l++){
        unsigned int vert0 = l;
        unsigned int vert1 = (l+1)%f.size();

        Vector3d tr;
        tr[0] = f[vert1][0] - f[vert0][0];
        tr[1] = f[vert1][1] - f[vert0][1];
        tr[2] = f[vert1][2] - f[vert0][2];

        Vector3d prod_t_t1 = t.cross(tr);
        Vector3d prod_t_p = (p-f[vert1]).cross(t);
        //Vector3d prod_t1_p = (f.vertices[vert1]-p).cross(tr);

        double alfa = (prod_t_p.dot(prod_t_t1))/(prod_t_t1.dot(prod_t_t1));
        //double beta = (prod_t1_p.dot(prod_t_t1))/(prod_t_t1.dot(prod_t_t1));

        if(alfa>=0-tol && alfa <=1+tol){
            Vector3d v;
            for(unsigned int i=0; i<3; i++){
                v[i] = f[vert0][i]*(alfa) + (1-alfa)*f[vert1][i];
            }
            intersectionsF.push_back(v);
        }

        if(intersectionsF.size()==2){
            if(abs(intersectionsF[0][0] - intersectionsF[1][0])<tol &&  abs(intersectionsF[0][1] - intersectionsF[1][1])<tol &&
                abs(intersectionsF[0][2] - intersectionsF[1][2])<tol ){
                intersectionsF.pop_back();
                //toRemove.push_back(l);
                continue;
            }
            //cout << "il primo poligono è tagliato" << endl;
            break;  //esco dal ciclo perchè di nuovo non è posisbile che una retta intersechi un poligono in più di due punti
        }           // potremmo anche mettere questo if fuori dall'else e levare il break prima , come ci piace di più
    }
    return intersectionsF;
}


VectorXd PALUSolver(const MatrixXd& a, const VectorXd& b){
    VectorXd solutionPALU = a.fullPivLu().solve(b);
    return solutionPALU;

}

bool orderLen(const Trace &a, const Trace &b){
    return a.len >= b.len;
}

bool compareFirstElement(const Vector3d& a, const Vector3d& b) {
    return (a[0] - b[0])<tol;
}


//calcolo distanza tra due punti
double dist(Vector3d v1, Vector3d v2){
    double sum= 0.0;
    for(unsigned int j =0; j<3; j++){
        double d = (v1[j] - v2[j])*(v1[j] - v2[j]);
        sum += d;
    }
    return sqrt(sum);
}



//controllo se un punto appartiene ad un segmento
bool onSegment(const Vector3d& p, const Vector3d& a, const Vector3d& b){
    double t = 0.0;
    unsigned int usata = 0;
    for(unsigned int i=0; i<3; i++){
        if(abs(b[i]-a[i])>tol){
            t = (p[i] - a[i])/(b[i]-a[i]);
            usata = i;
        }
    }
    //if((t<0.0)||(t>1.0))
    //    return false;
    vector<bool> check(2,false);
    unsigned int index = 0;
    for(unsigned int i=0; i<3; i++){
        if(i!=usata){
            double y = a[i] + t*(b[i]-a[i]);
            //cout << "p[i]: " << p[i] << " y: " << y << endl;
            if(abs(p[i] - y)<tol){
                check[index] = true;
            }
            index++;
        }
    }
    return check[0] && check[1];
}


void intersezioniSuRetta(bool& bole, vector<Vector3d>& trace, vector<Vector3d>& s1, vector<Vector3d>& s2){
    if(s1[0][0]  <= s2[0][0]){
        if(abs(s1[1][0]-s2[0][0])>tol && (s1[1][0]  <= s2[0][0])){
            return;
        }
        else if((s1[1][0] >= s2[0][0])  && (s1[1][0] <= s2[1][0]) ){
            trace[0] = s2[0];
            trace[1] = s1[1];
            bole = true;
        }
        else if(s1[1][0] >= s2[1][0]){
            trace[0] = s2[0];
            trace[1] = s2[1];
            bole = true;
        }

    }else if(s2[0][0]  <= s1[0][0]){
        if(abs(s2[1][0]-s1[0][0])>tol && (s2[1][0]  < s1[0][0])){
            return;
        }
        else if((s2[1][0] >=  s1[0][0]) && (s2[1][0] <= s1[1][0])){
            trace[0] = s1[0];
            trace[1] = s2[1];
            bole = true;
        }
        else if((s2[1][0] >= s1[1][0])){
            trace[0] = s1[0];
            trace[1] = s1[1];
            bole = true;
        }
    }
}




bool checkIsNew(const Cell1d& c2d, const Vector3d& point, const PolygonalMesh& pm, unsigned int& id){
    for(unsigned int c0d: c2d.extremes){
        Vector3d v = (pm.MapCell0D.at(c0d).coordinates);
        if(abs(v[0]-point(0))<tol && abs(v[1]-point(1))<tol && abs(v[2]-point(2))<tol){
            id = c0d;
            return false;
        }
    }

    return true;
}

unsigned int id0D = 0;
unsigned int id1D = 0;
unsigned int id2D = 0;

bool cuttingfractures(Cell2d& f, const Trace& t, PolygonalMesh& polyMesh, list<Cell2d>& next){
    vector<bool> angolo(2,false);
    unsigned int idSame;
    bool firstCell2d = true;
    unsigned int wheretoinsert;
    vector<Cell1d> forming;
    vector<Cell0d> forming0d;
    vector<unsigned int> idLatitagliati;
    unsigned int countcut = 0;
    Cell2d c2new1;
    Cell2d c2new2;
    c2new1.id = id2D++;
    c2new2.id = id2D++;
    unsigned int countline = 0;
    unsigned int countvertici = 0;
    bool beenFalse = false;
    vector<unsigned int> posVert;
    for(Cell1d& c1d: f.Cell2DEdges){
        c1d = polyMesh.MapCell1D.at(c1d.id);
        if(!c1d.old){
            unsigned int vert0 = countline;
            unsigned int vert1 = (countline+1)%(f.numVert);
            countline++;
            Vector3d intersection;
            //Vector3d p1 = pm.MapCell0D.at(c1d.extremes[0]).coordinates;
            cout << "pre try taglio" << endl;
            if(intersLato(t, f.Cell2DVertices[vert0], f.Cell2DVertices[vert1], angolo[countcut], intersection, polyMesh)){
                countcut++;
                if(checkIsNew(c1d, intersection, polyMesh, idSame)){
                    cout << "siamo qui" << endl;
                    f.old = true;
                    polyMesh.MapCell2D.at(f.id) = f;
                    idLatitagliati.push_back(c1d.id);
                    //c1d.old = true;
                    //polyMesh.MapCell1D.at(c1d.id) = c1d;
                    Cell0d newcell0d(id0D, intersection);

                    cout <<"e siamo anche qui" << endl;
                    vector<unsigned int> ext1 = {f.Cell2DVertices[vert0].id, id0D};
                    Cell1d newcell1d1(id1D, ext1);
                    id1D++;
                    vector<unsigned int> ext2 = {id0D++,f.Cell2DVertices[vert1].id};
                    Cell1d newcell1d2(id1D, ext2);
                    id1D++;
                    cout << "create celle 1d nuove" << endl;

                    c1d.tobecome.push_back(newcell1d1);
                    c1d.tobecome.push_back(newcell1d2);

                    /*newcell1d1.touched2D = c1d.touched2D;
                    newcell1d1.touched2D.push_back(c2new1.id); //havend't defined a buch  of this here yet
                    newcell1d2.touched2D = c1d.touched2D;
                    newcell1d2.touched2D.push_back(c2new2.id);*/
                    newcell0d.touched2D = c1d.touched2D;
                    newcell0d.touched2D.push_back(c2new1.id);
                    newcell0d.touched2D.push_back(c2new2.id);
                    forming0d.push_back(newcell0d);

                    if(firstCell2d){
                        newcell1d1.touched2D = c1d.touched2D;
                        newcell1d1.touched2D.push_back(c2new1.id); //havend't defined a buch  of this here yet
                        newcell1d2.touched2D = c1d.touched2D;
                        newcell1d2.touched2D.push_back(c2new2.id);

                        c2new1.Cell2DVertices.push_back(f.Cell2DVertices[vert0]); //quella da 200 si rompe qui
                        cout << "entrato" << endl;
                        c2new1.Cell2DVertices.push_back(newcell0d);
                        c2new1.Cell2DEdges.push_back(newcell1d1);

                        c2new2.Cell2DVertices.push_back(newcell0d);
                        c2new2.Cell2DEdges.push_back(newcell1d2);

                        wheretoinsert = c2new1.Cell2DEdges.size();
                        Cell1d tempCell;
                        c2new1.Cell2DEdges.push_back(tempCell);
                        firstCell2d = false;
                        beenFalse = true;

                        /*for(const unsigned int& c2did: c1d.touched2D){
                            if(!(polyMesh.MapCell2D.at(c2did)).old && c2did!=c2new1.id){
                                next.push_back(polyMesh.MapCell2D.at(c2did));
                            }
                        }*/
                    }
                    else{
                        newcell1d1.touched2D = c1d.touched2D;
                        newcell1d1.touched2D.push_back(c2new2.id); //havend't defined a buch  of this here yet
                        newcell1d2.touched2D = c1d.touched2D;
                        newcell1d2.touched2D.push_back(c2new1.id);
                        /*newcell0d.touched2D = c1d.touched2D;
                        newcell0d.touched2D.push_back(c2new2.id);
                        newcell0d.touched2D.push_back(c2new1.id);*/
                        c2new2.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
                        cout << "entrato nell'else"<<endl;
                        c2new2.Cell2DVertices.push_back(newcell0d);
                        c2new2.Cell2DEdges.push_back(newcell1d1);

                        c2new1.Cell2DVertices.push_back(newcell0d);
                        c2new1.Cell2DEdges.push_back(newcell1d2);
                        firstCell2d = true;

                        /*for(const unsigned int& c2did: c1d.touched2D){
                            if(!(polyMesh.MapCell2D.at(c2did)).old && c2did!=c2new2.id){
                                next.push_back(polyMesh.MapCell2D.at(c2did));
                            }
                        }*/
                    }
                    forming.push_back(newcell1d1);
                    forming.push_back(newcell1d2);

                    /*for(const unsigned int& c2did: c1d.touched2D){
                        if(!(polyMesh.MapCell2D.at(c2did)).old && c2did!=c2new1.id && c2did!=c2new2.id){
                            next.push_back(polyMesh.MapCell2D.at(c2did));
                        }
                    }*/
                    /*for(const unsigned int& c2did: polyMesh.MapCell0D.at(f.Cell2DVertices[vert0].id).touched2D){
                        if(!(polyMesh.MapCell2D.at(c2did)).old && c2did!=c2new1.id && c2did!=c2new2.id)
                            next.push_back(polyMesh.MapCell2D.at(c2did));
                    }
                    for(const unsigned int& c2did: polyMesh.MapCell0D.at(f.Cell2DVertices[vert1].id).touched2D){
                        if(!(polyMesh.MapCell2D.at(c2did)).old && c2did!=c2new1.id && c2did!=c2new2.id)
                            next.push_back(polyMesh.MapCell2D.at(c2did));
                    }*/
                }
                else{
                    if(firstCell2d){
                        posVert.push_back(idSame);
                        if(f.Cell2DVertices[vert0].id==idSame){
                            c2new1.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
                            c2new2.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
                            if(!beenFalse){
                                c2new2.Cell2DEdges.push_back(c1d);
                                firstCell2d = false;
                                beenFalse = true;

                            }else{
                                c2new1.Cell2DEdges.push_back(c1d);
                            }
                            /*wheretoinsert = c2new1.Cell2DEdges.size();
                            Cell1d tempCell;
                            c2new1.Cell2DEdges.push_back(tempCell);*/
                        }
                        else{
                            c2new1.Cell2DEdges.push_back(c1d);
                            c2new1.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
                            wheretoinsert = c2new1.Cell2DEdges.size();
                            Cell1d tempCell;
                            c2new1.Cell2DEdges.push_back(tempCell);
                        }
                        //firstCell2d = false;


                        /*for(const unsigned int& c2did: polyMesh.MapCell0D.at(f.Cell2DVertices[vert0].id).touched2D){
                            if(!(polyMesh.MapCell2D.at(c2did)).old)
                                next.push_back(polyMesh.MapCell2D.at(c2did));
                        }*/
                    }
                    else{
                        if(f.Cell2DVertices[vert0].id==idSame){
                            c2new1.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
                            c2new2.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
                            c2new2.Cell2DEdges.push_back(c1d);
                            /*wheretoinsert = c2new1.Cell2DEdges.size();
                            Cell1d tempCell;
                            c2new1.Cell2DEdges.push_back(tempCell);*/

                        }else{
                            c2new2.Cell2DEdges.push_back(c1d);
                            c2new2.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
                            firstCell2d = true;
                        }

                        /*for(const unsigned int& c2did: polyMesh.MapCell0D.at(f.Cell2DVertices[vert0].id).touched2D){
                            if(!(polyMesh.MapCell2D.at(c2did)).old)
                                next.push_back(polyMesh.MapCell2D.at(c2did));
                        }*/
                    }
                    /*for(const unsigned int& c2did: polyMesh.MapCell0D.at(f.Cell2DVertices[vert0].id).touched2D){
                        if(!(polyMesh.MapCell2D.at(c2did)).old && c2did!=c2new1.id && c2did!=c2new2.id)
                            next.push_back(polyMesh.MapCell2D.at(c2did));
                    }*/
                    /*Cell1d& c1pre = f.Cell2DEdges[countline-1];
                    for(const unsigned int& c2did: c1pre.touched2D){
                        if(!(polyMesh.MapCell2D.at(c2did)).old)
                            next.push_back(polyMesh.MapCell2D.at(c2did));
                    }*/
                }
            }
            else{
                countvertici++;
                if(firstCell2d){
                    c2new1.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
                    c2new1.Cell2DEdges.push_back(c1d);
                    c1d.touched2D.push_back(c2new1.id);
                    polyMesh.MapCell1D.at(c1d.id) = c1d;
                    polyMesh.MapCell0D.at(f.Cell2DVertices[vert0].id).touched2D.push_back(c2new1.id);

                }
                else{
                    c2new2.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
                    c2new2.Cell2DEdges.push_back(c1d);
                    c1d.touched2D.push_back(c2new2.id);
                    polyMesh.MapCell1D.at(c1d.id) = c1d;
                    polyMesh.MapCell0D.at(f.Cell2DVertices[vert0].id).touched2D.push_back(c2new2.id);
                }
            }
        }
    }
    if(!f.old)
        return false;
    if(countcut>1){   // tengo il conto delle intersezioni per evitare di creare degli oggetti nuovi qual'ora la traccia intersecasse la figure solo in un vertice
        Cell1d c1middle;
        c1middle.id = id1D++;
        if(posVert.size()==1){
            c1middle.extremes = {posVert[0], id0D-1};
        }
        else if(posVert.size()==2){
            c1middle.extremes = {posVert[0], posVert[1]};
        }
        else{
            c1middle.extremes={id0D-1, id0D-2};
        }
        c2new1.Cell2DEdges[wheretoinsert] = c1middle;
        c2new2.Cell2DEdges.push_back(c1middle);

        c1middle.touched2D.push_back(c2new1.id);
        c1middle.touched2D.push_back(c2new2.id);

        forming.push_back(c1middle);
        c2new1.numVert = c2new1.Cell2DVertices.size();
        c2new2.numVert = c2new2.Cell2DVertices.size();

        //newCell.push_back(c2new1);
        //newCell.push_back(c2new2);
        polyMesh.addingStuff(angolo, idLatitagliati, forming, forming0d, c2new1, c2new2, id2D, next);

        return true;
    }
    else{
        id2D = id2D-2;
        if(f.old){
            f.old = false;
            polyMesh.MapCell2D.at(f.id) = f;
            id0D = id0D-1;
            id1D = id1D-2;
        }

        return false;
    }
}

/*vector<unsigned int> removing;
        vector<unsigned int> coordt = intersezionipoligonorettaLATI(t.retta, t.coordTrace[1], f, removing);
        vector<Vector3d> coordtVere;

        if(f.tips.count(t.id)>0 && !f.tips.at(t.id)) //
            coordtVere = t.coordTrace;
        else
            coordtVere = intersezionipoligonoretta(t.retta, t.coordTrace[1], f);

        if(coordt.empty() || coordt.size()==1){  //la retta non interseca il poligono, aggiungo solo la fracture gia esistente
            //newfractures.resize(1);
            //newfractures[0] = f;
            return false;
        }

        list<unsigned int> vertOldFrac = f.idvertici;
        list<unsigned int> latiOldFrac = f.idlati;

        /*unsigned int numNewVertici = f.NumVertices+2;
        vector<bool> angolo(2, false);
        for(const auto& vertice: f.vertices){
            for(unsigned int a = 0; a<2; a++){
                if(vertice==coordtVere[a]){
                    numNewVertici--;
                    angolo[a] = true;
                    coordt[a]--;
                }
            }
        }*/

/*numNewVert = f.numVert + 2 - removing.size();
        Cell2d prima;
        Cell2d seconda;

        vector<Vector3d> newVertVec(numNewVertici);
        unsigned int j=0;
        unsigned int index = 0;
        //unsigned int uno = 1;
        auto it = vertOldFrac.begin();
        auto itLati = latiOldFrac.begin();
        //auto itlativecchi = f.idlati.begin();
        for(int i=0; i< numNewVert ; i++){
            if(i==coordt[index]){  // mi trovo nel punto in cui devo inserire il vertice nuovo
                for(const auto& r: removing){ // se il vertice nuovo corrisponde ad uno che esisteva già va solo avanti con le iteraizoni
                    if(r==i)
                        continue;
                }

                Cell0d newCell0d(id0D++, coordt[index]);
                Cell1d newCell1d(id1d++, )

                else{

                    polyMesh.Cell0DId.push_back(polyMesh.numCell0D);  //aggiungo il vertice al vettore di id e alla mappa dei vertivi
                    polyMesh.MapCell0D[polyMesh.numCell0D] = coordtVere[index];
                    vertOldFrac.insert(it, polyMesh.numCell0D);  //lo inserisco anche nella lista con tutti i vertici in ordne

                    if(index==0){
                        advance(it, -1);  // torno indietro di uno se no poi sforo (it itera sui vertici originali che quindi potrebbero essere due in meno di quelle che ho)
                        advance(itLati, -1);
                        unsigned int vecchiolatoid = *itLati;
                        vector<unsigned int> vecchiolato = polyMesh.MapCell1D[vecchiolatoid];
                        polyMesh.MapCell1D[polyMesh.numCell1D] = {vecchiolato[0], polyMesh.numCell0D};

                        *itLati = polyMesh.numCell1D++;
                        //latiOldFrac.remove(*itLati); //rimuovo il lato che ho tagliato
                        //latiOldFrac.insert(itLati, polyMesh.numCell1D++); //inserisco metà del lato tagliato
                        advance(itLati, 1);
                        polyMesh.MapCell1D[polyMesh.numCell1D] = {polyMesh.numCell0D++, vecchiolato[1]};
                        latiOldFrac.insert(itLati, polyMesh.numCell1D++); //inserisco l'altra metà del lato tagliato

                        //advance(itLati, 1);
                    }

                    if(index==1){
                        advance(itLati, -1);
                        unsigned int vecchiolatoid = *itLati;
                        vector<unsigned int> vecchiolato = polyMesh.MapCell1D[vecchiolatoid];
                        polyMesh.MapCell1D[polyMesh.numCell1D] = {vecchiolato[0], polyMesh.numCell0D};
                        advance(itLati, -1);
                        *itLati = polyMesh.numCell1D++;

                        //latiOldFrac.insert(itLati, polyMesh.numCell1D++);
                        advance(itLati, 1);
                        polyMesh.MapCell1D[polyMesh.numCell1D] = {polyMesh.numCell0D, polyMesh.numCell0D-1};
                        latiOldFrac.insert(itLati, polyMesh.numCell1D++);
                        //advance(itLati, 1);
                        polyMesh.MapCell1D[polyMesh.numCell1D] = {polyMesh.numCell0D++, vecchiolato[1]};
                        latiOldFrac.insert(itLati, polyMesh.numCell1D++);
                        //advance(itLati, -1);
                    }

                }
                //f.MapCell2DVertices.insert(i, polyMesh.numCell0D++);
                newVertVec[i] = coordtVere[index++];
            }
            else{
                newVertVec[i] = f.vertices[j++];
            }
            advance(it,1);
            advance(itLati, 1);
        }


        for(unsigned int i=0; i<coordt[0]; i++){
            unsigned int coda = vertOldFrac.front();
            vertOldFrac.pop_front();
            vertOldFrac.push_back(coda);

            unsigned int codalato = latiOldFrac.front();
            latiOldFrac.pop_front();
            latiOldFrac.push_back(codalato);
        }


        Fracture newFrac1, newFrac2;
        newFrac1.NumVertices = coordt[1] - coordt[0] + 1;
        newFrac1.vertices.resize(newFrac1.NumVertices);
        newFrac2.NumVertices = newVertVec.size() - newFrac1.NumVertices + 2;
        newFrac2.vertices.resize(newFrac2.NumVertices);
        unsigned int countlati2=0;


        auto iterVertici = vertOldFrac.begin();
        auto iterLati = latiOldFrac.begin();
        for(int pv = coordt[0]; pv<=coordt[1]; pv++){
            //cout << "punto prima figura" << endl;
            newFrac1.vertices[pv-coordt[0]] = newVertVec[pv];
            newFrac1.idvertici.push_back(*iterVertici);
            newFrac1.idlati.push_back(*iterLati);
            advance(iterLati, 1);

            if(pv!=coordt[1]){
                advance(iterVertici, 1);
            }
        }


        for( int pv = coordt[1]; pv < newVertVec.size(); pv++){
            //cout << "punto seconda figura" << endl;
            newFrac2.vertices[pv-coordt[1]] = newVertVec[pv];
            newFrac2.idvertici.push_back(*iterVertici);
            newFrac2.idlati.push_back(*iterLati);
            advance(iterVertici, 1);
            advance(iterLati, 1);
            countlati2++;
        }

        for(int pv = 0; pv <= coordt[0]; pv++){
            //.cout << "punto seconda figura" << endl;
            newFrac2.vertices[pv + (newVertVec.size() - coordt[1])] = newVertVec[pv];

            if(pv==coordt[0]){
                newFrac2.idvertici.push_back(vertOldFrac.front());
                advance(iterLati, -countlati2-1);
                newFrac2.idlati.push_back(*iterLati);
            }
            else{
                newFrac2.idvertici.push_back(*iterVertici);
                advance(iterVertici, 1);
                newFrac2.idlati.push_back((*iterLati));
                advance(iterLati, 1);
                countlati2++;
            }
        }
        //advance(iterLati, countlati2);
        //newFrac2.idlati.push_back(*iterLati);


        newfractures[0] = newFrac1;
        newfractures[1] = newFrac2;

        return newfractures;
    }*/


// importo i dati dai file
bool ImportFR_data(const string &filename, FractureMesh& mesh)
{
    //FractureMesh mesh;
    ifstream file;
    file.open(filename);
    if(file.fail())
        return false;
    string line;
    int N = 0;
    int i = 0;
    while (!file.eof())
    {
        getline(file, line);

        // Skip Comment Line
        if(line[0] != '#')
            break;
    }

    istringstream converter(line);
    converter.str(line);
    converter >> N; // stampo il numero di fratture nel file
    mesh.NumFractures = N;
    mesh.FractureId.resize(N);
    mesh.MapFractures.resize(N);

    char sep;
    unsigned int id = 0;
    unsigned int numvertices = 0;
    for (i=0; i < N; i++)
    {
        Fracture f;
        getline(file,line);
        getline(file, line); //leggo l'identificatore
        istringstream converter(line);
        converter >> id >> sep >> numvertices;
        f.id = id;
        mesh.FractureId[i]=id;
        f.NumVertices = numvertices;
        getline(file,line); // leggo riga vertici senza stamparla (vediamo se farlo tutto insieme le eliminazioni)

        Vector3d bary;
        f.vertices.resize(numvertices);

        for(unsigned int j=0; j<3; j++){
            getline(file, line);
            replace(line.begin(),line.end(),';',' ');
            istringstream converter(line);
            bary[j] = 0.0;

            for(unsigned int t =0; t<numvertices; t++){
                double a = 0.0;
                converter >> a;
                f.vertices[t][j] = a;
                bary[j]+=a;
            }
            bary[j] = bary[j]/numvertices;
        }

        f.barycentre = bary;
        mesh.MapFractures[id]=f;

    }

    file.close();
    return true;

}



void findIntersections(FractureMesh &mesh){
    for(unsigned int id = 0; id<mesh.NumFractures; id++){
        Fracture &f = mesh.MapFractures[id];
        Vector3d planeF = f.calcoloPiano();
        //double dF = -(planeF[0] * f.vertices[0][0]) - (planeF[1] * f.vertices[0][1]) - (planeF[2] * f.vertices[0][2]);
        //f.plane.reserve(4);

        //ciclo sui poigoni successivi
        for(unsigned int i=id+1; i<mesh.NumFractures; i++){
            bool inter = false;
            vector<Vector3d> trace = {};
            trace.resize(2);
            Vector3d t = {};
            Vector3d p = {};
            double distanzainF = 0;
            double distanzainFC = 0;
            Fracture &fConf = mesh.MapFractures[i];
            f.isOnEdge = false;
            fConf.isOnEdge =false;
            Vector3d planeFConf = fConf.calcoloPiano();
            //double dFConf = -(planeFConf[0] * fConf.vertices[0][0]) - (planeFConf[1] * fConf.vertices[0][1]) - (planeFConf[2] * fConf.vertices[0][2]);

            if((planeF.cross(planeFConf))[0]==0 && (planeF.cross(planeFConf))[1]==0 && (planeF.cross(planeFConf))[2]==0){
                if(abs(planeF[0]*fConf.plane[3] - planeFConf[0]*f.plane[3])<tol){  //sono complanari
                    for(unsigned int j=0; j<f.NumVertices; j++){
                        unsigned int vert0 = j;
                        unsigned int vert1 = (j+1)%f.NumVertices ;

                        Vector3d retta1 = {f.vertices[vert1][0]-f.vertices[vert0][0],f.vertices[vert1][1]-f.vertices[vert0][1],f.vertices[vert1][2]-f.vertices[vert0][2]};

                        for(unsigned int k=0; k<fConf.NumVertices; k++){
                            unsigned int vert2 = k;
                            unsigned int vert3 = (k+1)%fConf.NumVertices ;

                            vector<Vector3d> interFC;
                            if(sameLine(retta1, f.vertices[0], fConf.vertices, interFC)){
                                vector<Vector3d> seg1 = {f.vertices[vert0], f.vertices[vert1]};
                                vector<Vector3d> seg2 = {fConf.vertices[vert2], fConf.vertices[vert3]};
                                intersezioniSuRetta(inter, trace, seg1, seg2);
                                if(inter){
                                    f.isOnEdge = true;
                                    fConf.isOnEdge = true;
                                    distanzainF = dist(f.vertices[vert0], f.vertices[vert1]);
                                    distanzainFC = dist(fConf.vertices[vert2], fConf.vertices[vert3]);
                                }
                                break;
                            }
                        }
                    }
                }
                else{ //sono paralleli
                    continue;
                }
            }

            else if((dist(f.barycentre, fConf.barycentre) - (f.maxDist() + fConf.maxDist())) > tol ){
                // sicuramente non si intersecano perchè troppo lontane
                continue;

            }
            else { //se non sono complanari calcolo le intersezioni tra i piani
                t = planeF.cross(planeFConf);
                /*double nn = t.norm();
                for(unsigned int q=0; q<3; q++){
                    t[q] = t[q]/nn;
                }*/
                Matrix3d A;
                A.row(0) = planeF;
                A.row(1) = planeFConf;
                A.row(2) = t;
                Vector3d b;
                b[0] = -f.plane[3];
                b[1] = -fConf.plane[3];
                b[2] = 0.0;

                //cout << "b[0] == d di F " << b[0] << endl << "b[1] == d di fconf " << b[1] << endl;

                if(A.determinant() == 0){
                    break; //di nuovo il caso di complanarità, controlliamo solo per sicurezza
                }

                p = PALUSolver(A, b);

                //calcolo le interseioni dei lati della fracture con la retta di interseione dei piani
                vector<Vector3d> intersectionsF = intersezionipoligonoretta(t, p, f.vertices, f.isOnEdge);

                if(intersectionsF.size()==0){ // se nessuno dei lati del poligono interseca la retta allora sicuramente i due poligono in esame
                    // non si iintersecano -> passo a confromtarlo con un altro poligono
                    //cout << "Le figure " << id <<" e " << i <<" non si intersecano perche' la prima non e' attraversata dalla retta di intersezione tra i piani" << endl;
                    continue;
                }

                //calcolo le intersezioni della retta con l'altra fracture
                vector<Vector3d> intersectionsFC = intersezionipoligonoretta(t,p,fConf.vertices, fConf.isOnEdge);

                if(intersectionsFC.size()==0){
                    continue;
                }

                // calcolo le interseioni tra i due segmenti tovati
                sort(intersectionsF.begin(), intersectionsF.end(), compareFirstElement);
                sort(intersectionsFC.begin(), intersectionsFC.end(), compareFirstElement);

                intersezioniSuRetta(inter, trace, intersectionsF, intersectionsFC);

                if(inter){
                    distanzainF = dist(intersectionsF[0],intersectionsF[1]);
                    distanzainFC = dist(intersectionsFC[0], intersectionsFC[1]);
                }


            }

            if(inter){
                Trace newTrace;
                newTrace.coordTrace.resize(2);
                newTrace.fraId[0] =id;
                newTrace.fraId[1] = i;
                newTrace.coordTrace = trace;
                newTrace.retta = t;
                newTrace.p = p;
                defNewTrace(newTrace, distanzainF, distanzainFC, f, fConf, mesh);
                f.onEdge[newTrace.id] = f.isOnEdge;
                fConf.onEdge[newTrace.id] = fConf.isOnEdge;
            }
        }

        f.listPas.sort(orderLen);
        f.listNonpas.sort(orderLen);
    }
}

vector<PolygonalMesh> newpolygon(FractureMesh& mesh){
    vector<PolygonalMesh> newfigures;

    for(const Fracture& f: mesh.MapFractures){
        cout << "taglio frac:"  << f.id << endl;
        id0D=0;
        id1D=0;
        id2D=0;
        PolygonalMesh pm;

        Cell2d firstCell2d;
        firstCell2d.convertFracture(f, id0D, id1D, id2D);
        pm.addFirstCell2d(firstCell2d);

        list<Cell2d> next;
        //vector<Cell2d> newfrac;
        //newfrac.push_back(firstCell2d);
        for(const Trace& trace: f.listPas){
            unsigned int i =0;
            for(unsigned int& cc2d: pm.Cell2DId){
                Cell2d &cc = pm.MapCell2D.at(cc2d);
                if(!cc.old && cuttingfractures(cc, trace, pm, next)){
                    //newfrac.erase(i);
                    while(!next.empty()){
                        if(!next.front().old)
                        //Cell2d& n = *(next.front());
                            bool cutted = cuttingfractures(next.front(), trace, pm, next);
                            //newfrac.erase(next.front());
                        next.pop_front();
                    }
                    break;
                }

            }

        }

        for(const Trace& trace: f.listNonpas){
            vector<Vector3d> copiacordiTrace = trace.coordTrace;
            sort(copiacordiTrace.begin(), copiacordiTrace.end(), compareFirstElement);
            //cout << "trace (non passante) id:" << trace.id << endl;
            //> newfractures = {};
            for(unsigned int& cc2d: pm.Cell2DId){
                Cell2d &cc = pm.MapCell2D.at(cc2d);
                if(!cc.old){
                    /*Fracture ff;*/
                    vector<Vector3d> ver;
                    for(Cell0d& c0: cc.Cell2DVertices){
                        ver.push_back(c0.coordinates);
                    }
                    //ff.NumVertices = ff.vertices.size();
                    bool estremiTracciaInside = false;
                    bool onEdge = false;
                    vector<Vector3d> interRetta = intersezionipoligonoretta(trace.retta, trace.p, ver, onEdge);

                    if(interRetta.size()==2){
                        sort(interRetta.begin(), interRetta.end(), compareFirstElement);
                        vector<Vector3d> estremi(2);
                        intersezioniSuRetta(estremiTracciaInside, estremi, interRetta, copiacordiTrace);
                    }

                    if(estremiTracciaInside && !f.onEdge.at(trace.id)){
                        if(cuttingfractures(cc, trace, pm, next)){
                            //newfrac.remove(cc);
                            while(!next.empty()){
                                if(!next.front().old){
                                    vector<Vector3d> ver;
                                    for(Cell0d& c0: next.front().Cell2DVertices){
                                        ver.push_back(c0.coordinates);
                                    }
                                    //ff.NumVertices = ff.vertices.size();
                                    bool estremiTracciaInside2 = false;
                                    bool onEdge = false;
                                    vector<Vector3d> interRetta = intersezionipoligonoretta(trace.retta, trace.p, ver, onEdge);

                                    if(interRetta.size()==2){
                                        sort(interRetta.begin(), interRetta.end(), compareFirstElement);
                                        vector<Vector3d> estremi(2);
                                        intersezioniSuRetta(estremiTracciaInside2, estremi, interRetta, copiacordiTrace);
                                    }
                                    if(estremiTracciaInside2)
                                        bool cutted = cuttingfractures(next.front(), trace, pm, next);
                                }
                                next.pop_front();
                            }
                            break;
                        }

                    }
                    /*else{
                        Cell2d newcell;
                        newcell.id = id2D;
                        for(Cell1d& c1d: cc.Cell2DEdges){

                            if (c1d.old){
                                cout << "c1 old figura non tagliata" << endl;
                                newcell.Cell2DEdges.push_back(c1d.tobecome[0]);
                                newcell.Cell2DEdges.push_back(c1d.tobecome[1]);

                                newcell.Cell2DVertices.push_back(pm.MapCell0D.at(c1d.tobecome[0].extremes[0]));
                                newcell.Cell2DVertices.push_back(pm.MapCell0D.at(c1d.tobecome[0].extremes[1]));
                                newcell.Cell2DVertices.push_back(pm.MapCell0D.at(c1d.tobecome[1].extremes[1]));

                                pm.Cell2DId.push_back(id2D++);
                                cc.old = true;
                                cout << "tutto apporto" << endl;
                             }
                            /*else{
                                    newcell.Cell2DEdges.push_back(c1d);
                                    newcell.Cell2DVertices.push_back(pm.MapCell0D.at(c1d.extremes[0]));
                                    newcell.Cell2DVertices.push_back(pm.MapCell0D.at(c1d.extremes[1]));

                                }*/
                    //}
                    //}
                }
            }
        }
        newfigures.push_back(pm);
    }
    return newfigures;
}

void printingPolygonMesh(const vector<PolygonalMesh>& vpm, const string& file){
    string filepath = "newPolygons_" + file.substr(4,file.length()-4);
    ofstream outfile(filepath);
    if (outfile.is_open()) {
        for(const PolygonalMesh& pm: vpm){
            outfile << "New Polygonal Mesh" << endl;
            unsigned int countCell1d = 0;
            unsigned int countCell2d = 0;
            outfile << "Cell0Ds: " << endl << "Id; X; Y; Z;" << endl;
            for(const auto& d: pm.Cell0DId){
                outfile << d << "; " << pm.MapCell0D.at(d).coordinates[0] << "; " << pm.MapCell0D.at(d).coordinates[1] << "; " << pm.MapCell0D.at(d).coordinates[2] << endl;
            }
            outfile << "Cell0Ds: " << endl << "Id; IdVert1; IdVert2" << endl;
            for(const auto& e: pm.Cell1DId){
                if(!pm.MapCell1D.at(e).old){
                    outfile << e << "; " << pm.MapCell1D.at(e).extremes[0] << "; " << pm.MapCell1D.at(e).extremes[1] << endl;
                    countCell1d++;
                }
            }
            for(const auto &f: pm.Cell2DId){
                //if(!pm.MapCell2D.at(f).old){
                    countCell2d++;
                    outfile << "# Cell2dId; NumCell0d" << endl << f << "; " << pm.MapCell2DVertices.at(f).size() << endl << "# Cell0dId; X; Y; Z;" << endl;
                    for(const auto &t: pm.MapCell2DVertices.at(f)){
                        //auto it = pm.MapCell2DVertices.at(f).begin();
                        outfile << t << "; " ;
                        for(unsigned int k=0; k<3; k++)
                            outfile << scientific << setprecision(16)<< pm.MapCell0D.at(t).coordinates[k] << "; ";
                        outfile << endl;
                    }

                    outfile << "# Cell1ds; IdCell0ds; IdCell0ds;" << endl;
                    for(const auto &t: pm.MapCell2DEdges.at(f)){
                        //auto it = pm.MapCell2DEdges.at(f).begin();
                        outfile <<t << "; " ;
                        for(unsigned int k=0; k<2; k++)
                            outfile << scientific << setprecision(16)<< pm.MapCell1D.at(t).extremes[k] << "; ";
                        outfile << endl;
                    }
                //}
            }
            outfile << "Number of cell 0d: " << pm.Cell0DId.size() << endl;
            outfile << "Number of cell 1d: " << countCell1d << endl;
            outfile << "Number of cell 2d: " << countCell2d << endl;
            outfile << endl;

        }
        outfile.close();
    } else {
        cout << "Impossibile aprire il file." << endl;
    }
}


}
//Fracture f = coppia.second;

//unsigned int id0d = 0;
//unsigned int id1d = 0;
//unsigned int id2d = 0;
/*pm.Cell0DId.push_back(pm.numCell0D);
            pm.MapCell0D[pm.numCell0D] = f.vertices[0];
            f.idvertici.push_back(pm.numCell0D++);*/

/*for(unsigned int vert = 0; vert < f.NumVertices; vert ++){
                //unsigned int vert0 = vert-1;
                unsigned int vert1 = vert;
                pm.Cell0DId.push_back(pm.numCell0D);
                for(unsigned int x=0; x<3; x++){
                    pm.MapCell0D[pm.numCell0D][x] = f.vertices[vert1][x];
                }
                f.idvertici.push_back(pm.numCell0D++);

                pm.Cell1DId.push_back(pm.numCell1D);
                pm.MapCell1D[pm.numCell1D] = {pm.numCell0D-1, pm.numCell0D};
                f.idlati.push_back(pm.numCell1D++);

                if(pm.numCell0D==f.NumVertices)
                    pm.MapCell1D[pm.numCell1D] = {pm.numCell0D-1, 0};
                //pm.Cell1DId.push_bakc(pm.numCell1D);
                //pm.MapCell1D
            }

            /*map<unsigned int, Vector3d> MapCell0D = {};
            map<unsigned int, list<unsigned int>> MapCell1D = {};
            map<unsigned int, list<unsigned int>> MapCell2DEdges = {};*/

/*datagliare.push_back(firstCell2d);

            for(const Trace& trace: f.listPas){
                //cout << "la figura " << coppia.first << " ha numero tracce passanti: " << f.listPas.size() << endl;
                cout << "trace id:" << trace.id << endl;

                list<Fracture> newfractures = {};
                for(Cell2d& ff: datagliare){
                    if(!f.onEdge.at(trace.id)){
                        vector<Fracture> nfs = cuttingfractures(ff, trace, pm); //aggiorna la funzione (ti troan bool)
                                                                                // salvati tutte le cose nuove nella polyumesh
                                                                                // cicla du next dopo il primo taglio
                        //covrei poter fare la stessa cosa amnche per quelle non passanti
                        for(const Fracture& nf: nfs){
                            newfractures.push_back(nf);
                        }
                    }
                }
                datagliare.clear();
                datagliare = newfractures;

                /*unsigned int i=0;
                for (const Fracture& dd: datagliare){
                    cout << i++ << endl;
                    for (unsigned int h=0; h<dd.NumVertices; h++){
                        for(unsigned int t=0; t<3; t++){
                            cout << dd.vertices[h][t] << "  ";
                        }
                        cout << endl;
                    }
                    cout << endl;
                }*/

//}





/*for(const Fracture& nf: datagliare){
                cout << f.id << endl;
                pm.Cell2DId.push_back(id2d);
                pm.MapCell2DVertices[id2d] = nf.idvertici;
                pm.MapCell2DEdges[id2d] = nf.idlati;

                cout << id2d++ << ": " << endl;
                pm.numCell2D = id2d;
                //cout << nf.NumVertices<<endl;
                cout <<"celle 0d: " << endl;
                for(unsigned int j:  pm.MapCell2DVertices[id2d-1]){
                    cout << j << "  ";
                }
                cout << endl;

                cout << "celle 1d: " << endl;
                for(unsigned int k :  pm.MapCell2DEdges[id2d-1]){
                    cout << k << " ";
                }
                cout<<endl;
                /*if(abs(nf.vertices[j][k])<tol)
                            cout << scientific << setprecision(16) << 0.0 << " ";
                        else
                            cout << scientific << setprecision(16) << nf.vertices[j][k] << " ";
                    }
                    cout << endl;
                }
                id ++;*/
/*}
            datagliare.clear();
            newfigures[id++] = pm;
        }*/

















