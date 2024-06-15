/*#include "dfn.hpp"
#include "utils.hpp"
#include "poligonalMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
#include <algorithm>
#include <map>
#include <iomanip>
#include <list>

using namespace std;
using namespace Eigen;
using namespace Polygon;

namespace GeometryHelpers {

    double tol = 10 * numeric_limits<double>::epsilon();

bool  sameLine(const Vector3d& retta, const Fracture& f, vector<Vector3d> coordinate){
        double norma1 = retta.norm();
        Vector3d rettaN = {retta[0]/norma1, retta[1]/norma1, retta[2]/norma1};
        for(int k=0; k<f.NumVertices; k++){
            int vert0 = k;
            int vert1 = (k+1)%f.NumVertices;
            Vector3d retta2 = {f.vertices[vert1][0]-f.vertices[vert0][0],
                               f.vertices[vert1][1]-f.vertices[vert0][1],
                               f.vertices[vert1][2]-f.vertices[vert0][2]};
            double norma2 = retta2.norm();
            Vector3d retta2N = {retta2[0]/norma2, retta2[1]/norma2, retta2[2]/norma2};
            if(abs(abs(rettaN[0])-abs(retta2N[0]))<tol &&
                abs(abs(rettaN[1])-abs(retta2N[1]))<tol &&
                abs(abs(rettaN[1])-abs(retta2N[1]))<tol){
                coordinate.push_back(f.vertices[vert0]);
                coordinate.push_back(f.vertices[vert1]);
                return true;
            }
            else
                return false;
        }
    }




    bool intersLato(const Vector3d& t, Vector3d& p,  const Cell1d c1d, bool angolo, Vector3d& inters, PolygonalMesh& pm){
        Vector3d tr;
        Vector3d& p1 = pm.MapCell0D.at(c1d.extremes[0]);
        Vector3d& p2 = pm.MapCell0D.at(c1d.extremes[1]);
        tr[0] = p1[0] - p2[0];
        tr[1] = p1[1] - p2[1];
        tr[2] = p1[2] - p2[2];

        Vector3d prod_t_t1 = t.cross(tr);
        Vector3d prod_t_p = (p-p2).cross(t);
        //Vector3d prod_t1_p = (f.vertices[vert1]-p).cross(tr);

        double alfa = (prod_t_p.dot(prod_t_t1))/(prod_t_t1.dot(prod_t_t1));
        //double beta = (prod_t1_p.dot(prod_t_t1))/(prod_t_t1.dot(prod_t_t1));

        if(alfa>=0 && alfa <=1){
            for(unsigned int i=0; i<3; i++){
                inters[i] = p1[i]*alfa + (1-alfa)*p2[i];
            }
            if(abs(inters[0] - p1[0])<tol && abs(inters[1] - p1[1])<tol && abs(inters[2] - p1[2])<tol){
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

    /*vector<Vector3d> intersezionipoligonoretta(const Vector3d& t, const Vector3d p, Fracture& f){
        vector<Vector3d> intersectionsF = {};
        intersectionsF.reserve(2);

        //controllo che non la traccia non sia su un lato
        if(sameLine(t,f,intersectionsF)){
            f.isOnEdge = true;
            return intersectionsF;
        }
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
                Vector3d v;
                for(unsigned int i=0; i<3; i++){
                    v[i] = f.vertices[vert0][i]*alfa + (1-alfa)*f.vertices[vert1][i];
                }
                intersectionsF.push_back(v);
            }

            if(intersectionsF.size()==2){
                if(abs(intersectionsF[0][0] - intersectionsF[1][0])<tol &&  abs(intersectionsF[0][1] - intersectionsF[1][1])<tol &&
                    abs(intersectionsF[0][2] - intersectionsF[1][2])<tol ){
                    intersectionsF.pop_back();
                    toRemove.push_back(l);
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
        if((t<0.0)||(t>1.0))
            return false;
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

}*/
