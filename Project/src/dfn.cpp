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

namespace FractureLibrary {

double tol = 1000 * numeric_limits<double>::epsilon();
unsigned int idTrace = 0;

void defNewTrace(Trace& t, const double& d1, const double& d2, Fracture& f1, Fracture& f2, FractureMesh& fm){
    double lunghezzatraccia = dist(t.coordTrace[0],t.coordTrace[1]);
    t.len = lunghezzatraccia;
    t.id = idTrace;
    fm.MapTrace[idTrace++] = t;
    f1.numTrac++;
    f2.numTrac++;

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

bool  sameLine(const Vector3d& ret, const Vector3d& p, const vector<Vector3d>& f, vector<Vector3d>& coordinate){
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




bool intersLato(const Trace& tra, const Cell0d& c1, const Cell0d& c2, Vector3d& inters, const PolygonalMesh& pm){
    Vector3d t = tra.retta;
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
        return true;
    }
    else
        return false;
}



vector<Vector3d> intersezionipoligonoretta(const Vector3d& t, const Vector3d& p, vector<Vector3d>& f, bool& onEdge){
    vector<Vector3d> intersectionsF = {};
    intersectionsF.reserve(2);

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
            break;
        }
    }

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
    if(check[0] && check[1])
        return true;
    else
        return false;
}


void intersezioniSuRetta(bool& bole, vector<Vector3d>& trace, const vector<Vector3d>& s1, const vector<Vector3d>& s2){
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




bool checkIsNew(const vector<unsigned int>& c2d, const Vector3d& point, const PolygonalMesh& pm, unsigned int& id){
    for(unsigned int c0d: c2d){
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

void addNewVertAndEdg(bool& firstCell2d, unsigned int& wheretoinsert, Cell2d& c2new1, Cell2d& c2new2, bool& beenFalse,
                      vector<Cell1d>& forming,vector<Cell0d>& forming0d, Vector3d& intersection, Cell1d& c1d,
                      unsigned int& vert0, unsigned int& vert1, Cell2d& f){
    Cell0d newcell0d(id0D, intersection);
    //vector<unsigned int> ext1 = {f.Cell2DVertices[vert0].id, id0D};
    Cell1d newcell1d1(id1D, {f.Cell2DVertices[vert0].id, id0D});
    id1D++;
    //vector<unsigned int> ext2 = {id0D++,f.Cell2DVertices[vert1].id};
    Cell1d newcell1d2(id1D, {id0D++,f.Cell2DVertices[vert1].id});
    id1D++;

    newcell0d.touched2D = c1d.touched2D;
    newcell0d.touched2D.push_back(c2new1.id);
    newcell0d.touched2D.push_back(c2new2.id);

    if(firstCell2d){
        newcell1d1.touched2D = c1d.touched2D;
        newcell1d1.touched2D.push_back(c2new1.id);
        newcell1d2.touched2D = c1d.touched2D;
        newcell1d2.touched2D.push_back(c2new2.id);

        c2new1.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
        c2new1.Cell2DVertices.push_back(newcell0d);
        c2new1.Cell2DEdges.push_back(newcell1d1);

        c2new2.Cell2DVertices.push_back(newcell0d);
        c2new2.Cell2DEdges.push_back(newcell1d2);

        wheretoinsert = c2new1.Cell2DEdges.size();
        Cell1d tempCell;
        c2new1.Cell2DEdges.push_back(tempCell);
        firstCell2d = false;
        beenFalse = true;

    }
    else{
        newcell1d1.touched2D = c1d.touched2D;
        newcell1d1.touched2D.push_back(c2new2.id);
        newcell1d2.touched2D = c1d.touched2D;
        newcell1d2.touched2D.push_back(c2new1.id);

        c2new2.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
        c2new2.Cell2DVertices.push_back(newcell0d);
        c2new2.Cell2DEdges.push_back(newcell1d1);

        c2new1.Cell2DVertices.push_back(newcell0d);
        c2new1.Cell2DEdges.push_back(newcell1d2);
        firstCell2d = true;
    }
    forming0d.push_back(newcell0d);
    forming.push_back(newcell1d1);
    forming.push_back(newcell1d2);
    c1d.tobecome.push_back(newcell1d1);
    c1d.tobecome.push_back(newcell1d2);
}

void addNewEdg(bool& firstCell2d, unsigned int& wheretoinsert, Cell2d& c2new1, Cell2d& c2new2, bool& beenFalse, Cell1d& c1d,
                        Cell2d& f, PolygonalMesh& polyMesh, unsigned int& vert0, unsigned int& idSame){
    if(firstCell2d){
        c2new1.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
        c2new1.Cell2DVertices.push_back(polyMesh.MapCell0D.at(idSame));
        c2new2.Cell2DVertices.push_back(polyMesh.MapCell0D.at(idSame));
        c2new1.Cell2DEdges.push_back(polyMesh.MapCell1D.at(c1d.id).tobecome[1]);
        c2new2.Cell2DEdges.push_back(polyMesh.MapCell1D.at(c1d.id).tobecome[0]);

        wheretoinsert = c2new1.Cell2DEdges.size();
        Cell1d tempCell;
        c2new1.Cell2DEdges.push_back(tempCell);

        c1d.tobecome[1].touched2D.push_back(c2new1.id);
        c1d.tobecome[0].touched2D.push_back(c2new2.id);
        firstCell2d = false;
        beenFalse = true;

    }else{
        c2new2.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
        c2new2.Cell2DVertices.push_back(polyMesh.MapCell0D.at(idSame));
        c2new1.Cell2DVertices.push_back(polyMesh.MapCell0D.at(idSame));
        c2new2.Cell2DEdges.push_back(polyMesh.MapCell1D.at(c1d.id).tobecome[1]);
        c2new1.Cell2DEdges.push_back(polyMesh.MapCell1D.at(c1d.id).tobecome[0]);

        c1d.tobecome[1].touched2D.push_back(c2new2.id);
        c1d.tobecome[0].touched2D.push_back(c2new1.id);

        firstCell2d = true;
    }
    f.Cell2DVertices[vert0].touched2D.push_back(c2new1.id);
    f.Cell2DVertices[vert0].touched2D.push_back(c2new2.id);
}

void dividingExistingVert(const unsigned int& idSame, Cell2d& f, unsigned int& vert0, bool& firstCell2d, bool& beenFalse,
                        Cell2d& c2new1, Cell2d& c2new2, Cell1d& c1d, unsigned int& wheretoinsert ){
    if(firstCell2d){
        if(f.Cell2DVertices[vert0].id==idSame){
            c2new1.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
            c2new2.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);

            if(!beenFalse){
                c2new2.Cell2DEdges.push_back(c1d);
                c1d.touched2D.push_back(c2new2.id);
                firstCell2d = false;
                beenFalse = true;

            }else{
                c2new1.Cell2DEdges.push_back(c1d);
                c1d.touched2D.push_back(c2new1.id);
            }

        }
        else{
            c2new1.Cell2DEdges.push_back(c1d);
            c1d.touched2D.push_back(c2new1.id);
            c2new1.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
            if(!beenFalse){
                wheretoinsert = c2new1.Cell2DEdges.size();
                Cell1d tempCell;
                c2new1.Cell2DEdges.push_back(tempCell);
            }
        }
    }else{
        if(f.Cell2DVertices[vert0].id==idSame){
            c2new1.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
        }else{
            firstCell2d = true;
        }
        c2new2.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);

        c2new2.Cell2DEdges.push_back(c1d);
        c1d.touched2D.push_back(c2new2.id);
    }
}


bool cuttingfractures(Cell2d& f, const Trace& t, PolygonalMesh& polyMesh, list<Cell2d>& next){
    unsigned int idSame;
    bool firstCell2d = true;
    unsigned int wheretoinsert = 0;
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
    vector<unsigned int> idVertF;
    for(const Cell0d& p: f.Cell2DVertices){
        idVertF.push_back(p.id);
    }
    for(Cell1d& c1d: f.Cell2DEdges){
        c1d = polyMesh.MapCell1D.at(c1d.id);
        unsigned int vert0 = countline;
        unsigned int vert1 = (countline+1)%(f.numVert);

        Vector3d intersection;
        if(intersLato(t, f.Cell2DVertices[vert0], f.Cell2DVertices[vert1], intersection, polyMesh)){
            countcut++;
            f.old = true;
            if(checkIsNew(polyMesh.Cell0DId, intersection, polyMesh, idSame)){
                idLatitagliati.push_back(c1d.id);
                addNewVertAndEdg(firstCell2d, wheretoinsert, c2new1, c2new2, beenFalse, forming, forming0d, intersection, c1d, vert0, vert1, f);
            }
            else{
                posVert.push_back(idSame);
                unsigned int vertice;
                if(checkIsNew(idVertF, intersection, polyMesh, vertice)){
                    addNewEdg(firstCell2d, wheretoinsert, c2new1,c2new2,beenFalse,c1d,f, polyMesh, vert0, idSame);
                }else{
                    dividingExistingVert(idSame, f, vert0, firstCell2d, beenFalse, c2new1, c2new2, c1d, wheretoinsert);
                }
            }
        }
        else{
            countvertici++;
            if(firstCell2d){
                c2new1.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
                c1d.touched2D.push_back(c2new1.id);
                c2new1.Cell2DEdges.push_back(c1d);
                f.Cell2DVertices[vert0].touched2D.push_back(c2new1.id);
            }
            else{
                c2new2.Cell2DVertices.push_back(f.Cell2DVertices[vert0]);
                c1d.touched2D.push_back(c2new2.id);
                c2new2.Cell2DEdges.push_back(c1d);
                f.Cell2DVertices[vert0].touched2D.push_back(c2new2.id);
            }
        }
        countline++;

    }
    if(!f.old){
        id2D = id2D-2;
        return false;
    }
    if(countcut>1){   // tengo il conto delle intersezioni per evitare di creare degli oggetti nuovi qual'ora la traccia intersecasse la figure solo in un vertice
        unsigned int vert0 = 0;
        for(const Cell1d& cc1: f.Cell2DEdges){
            polyMesh.MapCell0D.at(f.Cell2DVertices[vert0].id) = f.Cell2DVertices[vert0];
            polyMesh.MapCell1D.at(cc1.id) = cc1;
            polyMesh.MapCell2D.at(f.id).Cell2DVertices[vert0] = f.Cell2DVertices[vert0];
            polyMesh.MapCell2D.at(f.id).Cell2DEdges[vert0] = cc1;
            if(cc1.old){
                polyMesh.MapCell1D.at(cc1.tobecome[0].id)  = cc1.tobecome[0];
                polyMesh.MapCell1D.at(cc1.tobecome[1].id)  = cc1.tobecome[1];
            }
            if(!cc1.tobecome.empty()){
                polyMesh.MapCell1D.at(cc1.id).old = true;
            }
            //polyMesh.MapCell1D.at(cc1.id) = cc1;
            vert0++;
        }
        //polyMesh.MapCell2D.at(f.id) = f;

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

        c1middle.touched2D.push_back(c2new1.id);
        c1middle.touched2D.push_back(c2new2.id);
        if(wheretoinsert==0)
            c2new1.Cell2DEdges.push_back(c1middle);
        else
            c2new1.Cell2DEdges[wheretoinsert] = c1middle;
        c2new2.Cell2DEdges.push_back(c1middle);

        forming.push_back(c1middle);
        c2new1.numVert = c2new1.Cell2DVertices.size();
        c2new2.numVert = c2new2.Cell2DVertices.size();

        polyMesh.MapCell2D.at(f.id) = f;
        polyMesh.addingStuff(idLatitagliati, forming, forming0d, c2new1, c2new2, id2D, next);

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

                Matrix3d A;
                A.row(0) = planeF;
                A.row(1) = planeFConf;
                A.row(2) = t;
                Vector3d b;
                b[0] = -f.plane[3];
                b[1] = -fConf.plane[3];
                b[2] = 0.0;

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

//funzione per controllare se la traccia non passante taglia una sotto-Cella2d
bool cuttedByNonPas(const vector<Vector3d>& copiacoordiTrace, const Cell2d& cc, const Fracture& f, const Trace& trace){

    vector<Vector3d> ver;
    for(const Cell0d& c0: cc.Cell2DVertices){
        ver.push_back(c0.coordinates);
    }
    bool estremiTracciaInside = false;
    bool onEdge = false;
    vector<Vector3d> interRetta = intersezionipoligonoretta(trace.retta, trace.p, ver, onEdge);

    if(interRetta.size()==2){
        sort(interRetta.begin(), interRetta.end(), compareFirstElement);
        vector<Vector3d> estremi(2);
        intersezioniSuRetta(estremiTracciaInside, estremi, interRetta, copiacoordiTrace);
    }

    if(estremiTracciaInside && !f.onEdge.at(trace.id))
        return true;
    else
        return false;
}

void splitOneEdg(unsigned int& id2D, Cell2d& toCut, PolygonalMesh& pm){
    Cell2d c2new3;
    c2new3.id = id2D;
    vector<unsigned int> listaidvertici;
    vector<unsigned int> listaidlati;
    toCut.old = true;
    unsigned int vert = 0;
    for(Cell1d& lati: toCut.Cell2DEdges){
        unsigned int idvert0 = toCut.Cell2DVertices[vert].id;
        if(pm.MapCell1D.at(lati.id).old){
            c2new3.Cell2DVertices.push_back(pm.MapCell0D.at(idvert0));
            listaidvertici.push_back(idvert0);
            Cell1d& c1 = pm.MapCell1D.at(lati.id).tobecome[0];
            Cell1d& c2 = pm.MapCell1D.at(lati.id).tobecome[1];
            c2new3.Cell2DVertices.push_back(pm.MapCell0D.at(c2.extremes[0]));
            c2new3.Cell2DEdges.push_back(c2);
            listaidlati.push_back(c2.id);
            c2new3.Cell2DEdges.push_back(c1);
            listaidlati.push_back(c1.id);
            listaidvertici.push_back(c2.extremes[0]);
            pm.MapCell1D.at(c1.id).touched2D.push_back(c2new3.id);
            pm.MapCell1D.at(c2.id).touched2D.push_back(c2new3.id);
            pm.MapCell1D.at(lati.id).tobecome[0].touched2D.push_back(c2new3.id);
            pm.MapCell1D.at(lati.id).tobecome[1].touched2D.push_back(c2new3.id);

        }else{
            c2new3.Cell2DVertices.push_back(pm.MapCell0D.at(idvert0));
            c2new3.Cell2DEdges.push_back(pm.MapCell1D.at(lati.id));
            listaidvertici.push_back(idvert0);
            listaidlati.push_back(lati.id);
            pm.MapCell1D.at(lati.id).touched2D.push_back(c2new3.id);
        }
        vert++;
    }

    pm.Cell2DId.push_back(id2D);
    c2new3.numVert = listaidvertici.size();
    pm.MapCell2D[id2D] = c2new3;
    pm.MapCell2DVertices[id2D] = listaidvertici;
    pm.MapCell2DEdges[id2D++] = listaidlati;
    pm.MapCell2D[toCut.id]=toCut;
}

vector<PolygonalMesh> newpolygon(FractureMesh& mesh){
    vector<PolygonalMesh> newfigures;

    for(const Fracture& f: mesh.MapFractures){
        id0D=0;
        id1D=0;
        id2D=0;
        PolygonalMesh pm;

        Cell2d firstCell2d;
        firstCell2d.convertFracture(f, id0D, id1D, id2D);
        pm.addFirstCell2d(firstCell2d);

        list<Cell2d> next;

        for(const Trace& trace: f.listPas){
            unsigned int i =0;
            for(unsigned int& cc2d: pm.Cell2DId){
                Cell2d &cc = pm.MapCell2D.at(cc2d);
                if(!cc.old && cuttingfractures(cc, trace, pm, next)){
                    while(!next.empty()){
                        if(!next.front().old)
                            bool cutted = cuttingfractures(next.front(), trace, pm, next);
                        next.pop_front();
                    }
                    break;
                }
            }
        }

        for(const Trace& trace: f.listNonpas){
            vector<Vector3d> copiacordiTrace = trace.coordTrace;
            sort(copiacordiTrace.begin(), copiacordiTrace.end(), compareFirstElement);

            for(unsigned int& cc2d: pm.Cell2DId){
                Cell2d &cc = pm.MapCell2D.at(cc2d);
                if(!cc.old){
                    if(cuttedByNonPas(copiacordiTrace, cc, f, trace)){
                        if(cuttingfractures(cc, trace, pm, next)){
                            while(!next.empty()){
                                Cell2d &toCut = next.front();
                                if(!toCut.old){
                                    if(cuttedByNonPas(copiacordiTrace, toCut, f, trace))
                                        bool cutted = cuttingfractures(toCut, trace, pm, next);
                                    else{
                                        splitOneEdg(id2D, toCut, pm);
                                    }
                                }
                                next.pop_front();
                            }
                            break;
                        }
                    }
                }
            }
        }
        newfigures.push_back(pm);
    }
    return newfigures;
}


//funzione per stampare tutte le polygonalMesh formate
void printingPolygonMesh(const vector<PolygonalMesh>& vpm, const string& file){
    string filepath = "newPolygons_" + file.substr(4,file.length()-4);
    ofstream outfile(filepath);
    if (outfile.is_open()) {
        for(const PolygonalMesh& pm: vpm){
            outfile << "New Polygonal Mesh" << endl;
            outfile << "Obtained from fracture " << pm.numFrac << endl;
            unsigned int countCell1d = 0;
            unsigned int countCell2d = 0;
            outfile << "Cell0Ds: " << endl << "Id; X; Y; Z;" << endl;
            for(const auto& d: pm.Cell0DId){
                outfile << d << "; " << pm.MapCell0D.at(d).coordinates[0] << "; " << pm.MapCell0D.at(d).coordinates[1] << "; " << pm.MapCell0D.at(d).coordinates[2] << endl;
            }
            outfile << "Cell1Ds: " << endl << "Id; IdVert1; IdVert2" << endl;
            for(const auto& e: pm.Cell1DId){
                if(!pm.MapCell1D.at(e).old){
                    outfile << e << "; " << pm.MapCell1D.at(e).extremes[0] << "; " << pm.MapCell1D.at(e).extremes[1] << endl;
                    countCell1d++;
                }
            }
            for(const auto &f: pm.Cell2DId){
                if(!pm.MapCell2D.at(f).old){
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
                }
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





