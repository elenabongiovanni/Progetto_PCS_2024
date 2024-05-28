#include "dfn.hpp"
#include "utils.hpp"
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

namespace FractureLibrary {

    double tol = 10 * numeric_limits<double>::epsilon();

    /*double findMin(const vector<double> v){
        auto minIter = min_element(v.begin(), v.end());

        if (minIter != v.end()) {
            double minValue = *minIter;
            return minValue;
        }  else {
            return -1;
        }
           /* std::cout << "Il valore minimo nel vettore è: " << minValue << std::endl;
        } else {
            std::cout << "Il vettore è vuoto." << std::endl;
        }
    }

    double findMax(const vector<double> v){
        auto maxElementIterator = max_element(v.begin(), v.end());

        if (maxElementIterator != v.end()) {
            double maxValue = *maxElementIterator;
            return maxValue;
            //std::cout << "Il valore massimo nel vettore è: " << maxValue << std::endl;
        } else {
            //std::cout << "Il vettore è vuoto." << std::endl;
            return -1;
        }
    }*/


    bool isPointInsideFracture(const Vector3d& p, const Fracture& f) {

        int n = f.vertices.size();
        /*vector<double> x(2);
        vector<double> y(2);
        vector<double> z(2);*/
        list<double> x;
        list<double> y;
        list<double> z;
        for(unsigned int i=0; i<n; i++){
            /*x[i] = f.vertices[i][0];
            y[i] = f.vertices[i][1];
            z[i] = f.vertices[i][2];*/
            x.push_back(f.vertices[i][0]);
            y.push_back(f.vertices[i][1]);
            z.push_back(f.vertices[i][2]);
        }
        x.sort();
        y.sort();
        z.sort();
        bool inside = false;

        if((p[0]-x.front())>=-tol && (p[0]- x.back())<=tol &&
            (p[1]- y.front())>=-tol && (p[1]-y.back())<=tol &&
            (p[2]-z.front())>=-tol && (p[2] - z.back())<=tol){
            inside = true;
        }

        /*for (int i = 0, j = n - 1; i < n; j = i++) {
            if ((f.vertices[i][1] > p[1]) != (f.vertices[j][1] > p[1]) &&
                p[0] < (f.vertices[j][0] - f.vertices[i][0]) * (p[1] - f.vertices[i][1]) / (f.vertices[j][1] - f.vertices[i][1]) + f.vertices[i][0] &&
                (p[2] < (f.vertices[j][2] - f.vertices[i][2]) * (p[1] - f.vertices[i][1]) / (f.vertices[j][1] - f.vertices[i][1]) + f.vertices[i][2])){
                inside = true;
            }
        }*/

        return inside;
    }

    Vector3d calcoloPiano(const Fracture& f){
        Vector3d lato1F = {};
        Vector3d lato2F = {};
        //calcolo 2 lati
        for(unsigned int j =0; j<3; j++){
            lato1F[j] = f.vertices[1][j] - f.vertices[0][j];
            lato2F[j] = f.vertices[3][j] - f.vertices[0][j];
        }
        //normalizzo i lati
        double lato1Fn = lato1F.squaredNorm();
        double lato2Fn = lato2F.squaredNorm();
        //calcolo il vettore normale e lo normalizzo
        Vector3d planeF = {};
        planeF[0] = (lato1F[1]*lato2F[2] - lato1F[2]*lato2F[1])/(lato1Fn*lato2Fn);
        planeF[1] = (lato1F[2]*lato2F[0] - lato1F[0]*lato2F[2])/(lato1Fn*lato2Fn);
        planeF[2] = (lato1F[0]*lato2F[1] - lato1F[1]*lato2F[0])/(lato1Fn*lato2Fn);
        return planeF;
    }

    vector<Vector3d> intersezionipoligonoretta(const Vector3d& t, const Vector3d p, const Fracture& f){
        vector<Vector3d> intersectionsF;
        intersectionsF.reserve(2);
        for(unsigned int l = 0; l < f.NumVertices; l++){
            unsigned int &vert0 = l;
            unsigned int vert1 = l+1;

            if(l==f.NumVertices-1){
                vert1 = l-(f.NumVertices-1);
            }

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
                //cout << "il primo poligono è tagliato" << endl;
                break;  //esco dal ciclo perchè di nuovo non è posisbile che una retta intersechi un poligono in più di due punti
            }           // potremmo anche mettere questo if fuori dall'else e levare il break prima , come ci piace di più
        }
        return intersectionsF;
    }

    vector<Fracture> cuttingfractures(const Fracture& f, const Trace& t){
        vector<Fracture> newfractures(2);

        vector<Vector3d> coordt = intersezionipoligonoretta(t.retta, t.coordTrace[1], f);
        if(coordt.empty()){
            newfractures.resize(1);
            newfractures.push_back(f);
            return newfractures;
        }

        //cout << coordt[0][0] << " " << coordt[0][1]<< " " << coordt[0][2] << endl;
        //cout << coordt[1][0] << " " << coordt[1][1]<< " " << coordt[1][2] << endl;

        list<Vector3d> newVert;
        //newVert.reserve(f.NumVertices+2);
        for(const Vector3d& cv: f.vertices){
            newVert.push_back(cv);
        }
        Vector2i posPunti= {-1, -1};
        unsigned int index = 0;
        for(Vector3d &punto: coordt){
            //cout << punto[0] << " " << punto[1] << " " << punto[2] << endl;

            for(unsigned int l = 0; l<f.NumVertices; l++){

                unsigned int vert0 = l;
                unsigned int vert1 = l+1;
                if(l==f.NumVertices-1){
                    vert1 = l-(f.NumVertices-1);
                }
                /*for(unsigned int i=0; i<3; i++){
                    if(abs(punto[i])<tol)
                        punto[i]=0;
                }*/

                if(onSegment(punto, f.vertices[vert0], f.vertices[vert1])){
                    //cout << "onsegment" << endl;
                    if(f.vertices[vert0][0]==punto[0] && f.vertices[vert0][1]==punto[1] && f.vertices[vert0][2]==punto[2])
                        posPunti[index] = vert0;
                    else if(f.vertices[vert1][0]==punto[0] && f.vertices[vert1][1]==punto[1] && f.vertices[vert1][2]==punto[2])
                        posPunti[index] = vert1;
                    else{
                        //cout << "siamo nel lato " << vert0 << "-" << vert1 << endl;
                        auto it = newVert.begin();
                        advance(it, vert0 + 1 + index); // Correggi l'inserimento
                        newVert.insert(it, punto);
                        posPunti[index] = vert0 + 1 + index;
                    }
                break;
                }
            }
            index++;
        }


        Fracture newFrac1, newFrac2;
        newFrac1.NumVertices = posPunti[1] - posPunti[0] + 1;
        newFrac1.vertices.resize(newFrac1.NumVertices);
        newFrac2.NumVertices = newVert.size() - newFrac1.NumVertices + 2;
        newFrac2.vertices.resize(newFrac2.NumVertices);

        vector<Vector3d> newVertVec(newVert.size());
        unsigned int i = 0;
        for(const Vector3d& v3 : newVert) {
            newVertVec[i]=v3;
            i++;
            //cout << "Adding vertex: " << v3[0] << " " << v3[1] << " " << v3[2] << endl;

        }
        //cout << "newVertVec size: " << newVertVec.size() << endl;
        //cout << "posPunti: " << posPunti[0] << ", " << posPunti[1] << endl;


        for(int pv = posPunti[0]; pv<=posPunti[1]; pv++){
            //cout << "punto prima figura" << endl;
            newFrac1.vertices[pv-posPunti[0]] = newVertVec[pv];
        }


        for( int pv = posPunti[1]; pv < newVert.size(); pv++){
            //cout << "punto seconda figura" << endl;
            newFrac2.vertices[pv-posPunti[1]] = newVertVec[pv];
        }

        for( int pv = 0; pv <= posPunti[0]; pv++){
            //.cout << "punto seconda figura" << endl;
            newFrac2.vertices[pv + (newVert.size() - posPunti[1])] = newVertVec[pv];
        }

        newfractures[0] = newFrac1;
        newfractures[1] = newFrac2;

        /*cout << "frac1:" << endl;
        for(int p = 0; p<newFrac1.NumVertices; p++){
            for(int k=0; k<3; k++){
                cout << newFrac1.vertices[p][k] << " ";
            }
            cout << endl;
        }

        cout << "frac2:" << endl;
        for(int p = 0; p<newFrac2.NumVertices; p++){
            for(int k=0; k<3; k++){
                cout << newFrac2.vertices[p][k] << " ";
            }
            cout << endl;
        }*/

        return newfractures;
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

    //calcolo distanza massima dal baricentro
    double maxDist(Fracture f){
        double r = 0.0;
        for(unsigned int i=0; i< f.NumVertices; i++){
            double d = dist(f.barycentre, f.vertices[i]);
            if(d>r)
                r = d;
        }
        return r;
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



    // importo i dati dai file
    bool ImportFR_data(const string &filename,
                       FractureMesh& mesh)
    {
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

        //cout << "numero di fratture: " << mesh.NumFractures << endl;

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

            //cout << "la fracture  " << mesh.FractureId[i] << " ha " << numvertices << " vertici" << endl;
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

    unsigned int idTrace = 0;

    void findIntersections(FractureMesh &mesh){
        for(unsigned int id = 0; id<mesh.NumFractures; id++){
        Fracture &f = mesh.MapFractures.at(id);

        //mesh.vecTrace.reserve(mesh.NumFractures*(mesh.NumFractures-1));
        Vector3d planeF = calcoloPiano(f);
        double dF = -(planeF[0] * f.vertices[0][0]) - (planeF[1] * f.vertices[0][1]) - (planeF[2] * f.vertices[0][2]);

        //ciclo sui poigoni successivi
        for(unsigned int i=id+1; i<mesh.NumFractures; i++){
            bool inter = false;
            vector<Vector3d> trace(2);
            //trace.resize(2);
            Vector3d t = {};
            double distanzainF = 0;
            double distanzainFC = 0;

            Fracture &fConf = mesh.MapFractures.at(i);

            Vector3d planeFConf = calcoloPiano(fConf);
            double dFConf = -(planeFConf[0] * fConf.vertices[0][0]) - (planeFConf[1] * fConf.vertices[0][1]) - (planeFConf[2] * fConf.vertices[0][2]);

            //cout <<"Il piano "<< id << " ha equazione "<< planeF[0] <<"x + " <<planeF[1] <<"y + " <<planeF[2] <<"z + " << dF << " = 0 "<<endl;
            //cout <<"Il piano "<< i << " ha equazione "<< planeFConf[0] <<"x + " <<planeFConf[1] <<"y + " <<planeFConf[2] <<"z + " << dFConf << "=0 "<<endl;

            if(planeF[0]==planeFConf[0] && planeF[1]==planeFConf[1] && planeF[2]==planeFConf[2] && dF == dFConf){ //condizone di complanarità
                //cout << "le figure "<< id <<" e " << i <<" sono complanari" << endl;
                /*for(int j=0; j<f.NumVertices; j++){
                    int &vert0 = j;
                    int vert1 = j+1;
                    if(j==f.NumVertices-1){
                        vert1 = j-f.NumVertices-1;
                    }
                    Vector3d retta(f.vertices[vert1][0]-f.vertices[vert0][0],f.vertices[vert1][1]-f.vertices[vert0][1],f.vertices[vert1][2]-f.vertices[vert0][2]);
                    for(int k=0; k<fConf.NumVertices; k++){
                        int &vert0 = k;
                        int vert1 = k+1;
                        if(k==fConf.NumVertices-1){
                            vert1 = k-3;
                        }

                }*/


                //continue;
            }

            else if(planeF[0]==planeFConf[0] && planeF[1]==planeFConf[1] && planeF[2]==planeFConf[2] && dF != dFConf){
                //cout << "le figure "<< id <<" e " << i <<" giacciono su piani paralleli" << endl;
                continue;
            }

            else if((dist(f.barycentre, fConf.barycentre) - (maxDist(f) + maxDist(fConf))) > tol ){
                //cout << "le figure "<< id <<" e " << i <<" sicuramente non si intersecano" << endl;
                continue;

            }
            else { //se non sono complanari calcolo le intersezioni tra i piani
                t = planeF.cross(planeFConf);

                Matrix3d A;
                A.row(0) = planeF;
                A.row(1) = planeFConf;
                A.row(2) = t;
                Vector3d b;
                b[0] = -dF;
                b[1] = -dFConf;
                b[2] = 0.0;

                /*if(A.determinant() == 0){
                    //cout << "Le figure "<< id <<" e " << i <<" il deteminante è nullo ora non mi ricordo il perchè" << endl;
                    break; //questo non è di nuovo il caso di complanarità???
                }*/

                Vector3d p = PALUSolver(A, b);

                //calcolo le interseioni dei lati della fracture con la retta di interseione dei piani
                vector<Vector3d> intersectionsF = intersezionipoligonoretta(t, p, f);

                if(intersectionsF.size()==0){ // se nessuno dei lati del poligono interseca la retta allora sicuramente i due poligono in esame
                                              // non si iintersecano -> passo a confromtarlo con un altro poligono
                    //cout << "Le figure " << id <<" e " << i <<" non si intersecano perche' la prima non e' attraversata dalla retta di intersezione tra i piani" << endl;
                    continue;
                }


                //calcolo le intersezioni della retta con l'altra fracture
                vector<Vector3d> intersectionsFC = intersezionipoligonoretta(t,p,fConf);

                if(intersectionsFC.size()==0){
                    //cout << "Le figure " << id <<" e " << i << " non si intersecano perche' la seconda non e' attraversata dalla retta di intersezione tra i piani" << endl;
                    continue;
                }

                // calcolo le interseioni tra i due segmenti tovati
                sort(intersectionsF.begin(), intersectionsF.end(), compareFirstElement);
                sort(intersectionsFC.begin(), intersectionsFC.end(), compareFirstElement);

                if(intersectionsF[0][0] <= intersectionsFC[0][0]){
                    if(abs(intersectionsF[1][0]-intersectionsFC[0][0])>tol && intersectionsF[1][0]<intersectionsFC[0][0]){
                        //cout <<"Le figure "<< id <<" e " << i << " non si intersecano" << endl;
                        continue;
                    }
                    else if(intersectionsF[1][0] >= intersectionsFC[0][0] && intersectionsF[1][0]<= intersectionsFC[1][0]){
                        trace[0]=intersectionsFC[0];
                        trace[1]=intersectionsF[1];
                        //cout << "LE DUE FIGURE " << id <<" e " << i << " SI INTERSECANO" << endl;
                        inter = true;

                    }
                    else if(intersectionsF[1][0]>=intersectionsFC[1][0]){                        
                        trace[0]=intersectionsFC[0];
                        trace[1]=intersectionsFC[1];
                        //cout << "LE DUE FIGURE " << id <<" e " << i <<" SI INTERSECANO" << endl;
                        inter = true;
                    }

                }else if(intersectionsFC[0][0]<=intersectionsF[0][0]){
                    if(abs(intersectionsFC[1][0]-intersectionsF[0][0])>tol && intersectionsFC[1][0]<intersectionsF[0][0]){
                        //cout <<"Le figure "<< id <<" e " << i << " non si intersecano" << endl;
                        continue;
                    }
                    else if(intersectionsFC[1][0]>=intersectionsF[0][0] && intersectionsFC[1][0]<=intersectionsF[1][0]){
                        trace[0]=(intersectionsF[0]);
                        trace[1]=(intersectionsFC[1]);
                        //cout << "LE DUE FIGURE " << id <<" e " << i << " SI INTERSECANO" << endl;
                        inter = true;

                    }
                    else if(intersectionsFC[1][0]>=intersectionsF[1][0]){
                        trace[0]=(intersectionsF[0]);
                        trace[1]=(intersectionsF[1]);
                        //cout << "LE DUE FIGURE " << id <<" e " << i << " SI INTERSECANO" << endl;
                        inter = true;
                    }
                }
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
                double lunghezzatraccia = dist(trace[0],trace[1]);
                newTrace.len = lunghezzatraccia;
                newTrace.id = idTrace;
                mesh.MapTrace[idTrace] = newTrace;
                idTrace++;
                f.numFrac++;
                fConf.numFrac++;
                newTrace.retta = t;

                if(abs(lunghezzatraccia-distanzainF)<tol){
                    //cout << "figura " << id << " e' passante"<< endl;
                    f.listPas.push_back(newTrace);
                    f.tips[newTrace.id] = false;

                }else{
                    f.listNonpas.push_back(newTrace);
                    f.tips[newTrace.id] = true;
                }

                if(abs(lunghezzatraccia-distanzainFC)<tol){
                    //cout << "figura " << i << " e' passante"<<endl;
                    fConf.listPas.push_back(newTrace);
                    fConf.tips[newTrace.id] = false;

                }else{
                    fConf.listNonpas.push_back(newTrace);
                    fConf.tips[newTrace.id] = true;
                }
            }
        }

        f.listPas.sort(orderLen);
        f.listNonpas.sort(orderLen);
        }
    }

    void printingtraces(FractureMesh& mesh, const string& file){
        string filepath = "traces_" + file.substr(4,file.length()-4);
        ofstream outfile(filepath);
        if (outfile.is_open()) {
            outfile << "#N Traces" << endl << mesh.MapTrace.size() << endl;
            outfile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
            for(const auto &t: mesh.MapTrace){
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

    void printingfractures(FractureMesh& mesh, const string& file){
        string filepath = "fractures_" + file.substr(4,file.length()-4);
        ofstream outfile(filepath);
        if (outfile.is_open()) {
            for(const auto &f: mesh.MapFractures){
                outfile << "# FractureId; NumTraces" << endl << f.id << "; " << f.numFrac << endl << "# TraceId; Tips; Lenght" << endl;
                for(const auto &t: f.tips){
                    outfile << t.first << "; ";
                    if(t.second){
                        outfile << "true; ";
                    }
                    else{
                        outfile << "false; ";
                    }
                    outfile << scientific << setprecision(16)<< mesh.MapTrace[t.first].len << endl;
                }
                outfile << endl;
            }
            outfile.close();
        } else {
            cout << "Impossibile aprire il file." << endl;
        }
    }

    /*PolygonalMesh newpolygon(FractureMesh& mesh){

        PolygonalMesh pm;
        unsigned int id = 0;
        //map<unsigned int, Fracture> MapNewFractures = {};

        for(const Fracture& f: mesh.MapFractures){
            //Fracture &f = mesh.MapFractures[chiave];

            cout << "taglio della frattura " << f.id << endl;

            list<Fracture> datagliare = {f};

            for(const Trace& trace: f.listPas){
                unsigned int lun = datagliare.size();

                for(unsigned int n=0; n<lun; n++){
                    Fracture &ff = datagliare.front();
                    vector<Fracture> nfs = cuttingfractures(ff, trace);
                    datagliare.pop_front();
                    cout << "dimensioni nfs: " << nfs.size() << endl;
                    //cout << "dimensioni new frac: " << newfractures.size() << endl;
                    for(const Fracture& nf: nfs){
                        datagliare.push_back(nf);

                    }

                    cout << "datagliare: " << datagliare.size() << endl;
                }



            }
            cout << "datagliare " << datagliare.size() << endl;
            cout << "inizio tracce non passanto" << endl;

            for(const Trace& trace: f.listNonpas){

                unsigned int lun = datagliare.size();

                for(unsigned int n=0; n<lun; n++){
                    Fracture &ff = datagliare.front();
                    //cout << ff.vertices[0][0] << " " << ff.vertices[0][1] << " " << ff.vertices[0][2] << endl << ff.vertices[1][0] << " " << ff.vertices[1][1] << " " << ff.vertices[1][2] << endl;
                    bool cutted = false;

                    for(unsigned int i=0; i<ff.NumVertices; i++){
                        unsigned int &vert0 = i;
                        unsigned int vert1 = i+1;
                        if(i==f.NumVertices-1){
                            vert1 = i-(f.NumVertices-1);
                        }
                        if(onSegment(trace.coordTrace[0], ff.vertices[vert0], ff.vertices[vert1]) || onSegment(trace.coordTrace[1], ff.vertices[vert0], ff.vertices[vert1])){
                            cutted = true;

                            cout << "datagliare: " << datagliare.size() << endl;
                            break;
                        }

                    }
                    if(cutted){
                        vector<Fracture> nfs = cuttingfractures(ff, trace);

                        for(const Fracture& nf: nfs){
                            datagliare.push_back(nf);
                        }
                        //cout << "added to the end: " << ff.vertices[0][0] << " " << ff.vertices[0][1] << " " << ff.vertices[0][2] << endl;
                        //continue;
                    }
                    else{
                        datagliare.push_back(ff);
                    }
                    datagliare.pop_front();

                }

            }

            unsigned int lun = datagliare.size();
            for(unsigned int n=0; n<lun; n++){
                Fracture &nf = datagliare.front();


                cout << id << ": " << endl;
                for(unsigned int j=0; j<nf.NumVertices; j++){
                    for(unsigned int k = 0; k<3; k++){
                        cout << nf.vertices[j][k] << " ";
                    }
                    cout << endl;
                }
                //datagliare.pop_front();
                id++;
            }
        }
        // dalle fratture ricavati lati e vertici -> trova modo per non mettere gloi stessi due volte

        return pm;
    }*/

    /*PolygonalMesh newpolygon(FractureMesh& mesh){

        PolygonalMesh pm;

        unsigned int id = 0;
        //map<unsigned int, Fracture> MapNewFractures = {};

        for(const Fracture& f: mesh.MapFractures){
            //Fracture &f = mesh.MapFractures[chiave];

            cout << "taglio della frattura " << f.id << endl;
            //vector<Fracture> &nfs = cuttingfractures(f, trace ,pm);
            /*map<unsigned int, Vector3d> MapCell0D = {};
            map<unsigned int, list<unsigned int>> MapCell1D = {};
            map<unsigned int, list<unsigned int>> MapCell2DEdges = {};*/
            //vector<Fracture> newfracture = {f};
            //list<Fracture> datagliare = {f};
            /*datagliare.resize(1);
            datagliare[0] = f;*/


            //for(const Trace& trace: f.listPas){
                //unsigned int newsize = 0;
                //cout << trace.id << endl;
                //cout << "la figura " << coppia.first << " ha numero tracce passanti: " << f.listPas.size() << endl;
                //vector<vector<Fracture>> vecCut(f.listPas.size());
                //newfractures.reserve(2*datagliare.size());
                //list<Fracture> newfractures;
               // unsigned int i = 0;
                /*for(const Fracture &ff: datagliare){
                    vector<Fracture> nfs = cuttingfractures(ff, trace);
                    cout << "dim nfs: " << nfs.size() << endl;
                    for(const Fracture& nf: nfs){
                        newfractures.push_back(nf);

                    }

                    //vecCut[i++]=nfs;
                    //newsize+=nfs.size();
                }
                datagliare = newfractures;
                /*unsigned int c = 0;
                unsigned int s = newfractures.size();
                for(unsigned int i=0;  i<s; i++){
                    if(newfractures[i] == Fracture{})
                        c++;
                }
                newfractures.resize(s-c);*/

                //datagliare.clear();
                //datagliare.resize(newfractures.size());
                //for(unsigned int i=0; i<newfractures.size();i++)
                //    datagliare[i]=newfractures[i];
                //newfractures.clear();
                //}
                //datagliare.resize(newfractures.size());
                /*datagliare.resize(newsize);
                unsigned int j = 0 ;
                for(const auto& cut:vecCut){
                    for(const Fracture& fra: cut){
                        datagliare[j++] = fra;
                    }
                }*/



            //}


            /*for(const Trace& trace: f.listNonpas){
                bool cutted = false;
                //unsigned int newsize = 0;
                //newfractures.clear();
                //cout << trace.id << endl;
                //cout << "la figura " << coppia.first << " ha numero tracce non passanti: " << f.listNonpas.size() << endl;
                //vector<Fracture> newfractures = {};
                //newfractures.reserve(datagliare.size()*2);
                list<Fracture> newfractures;
                //vector<vector<Fracture>> vecCut(f.listPas.size());
                //newfractures.reserve(2*datagliare.size());
                //unsigned int k = 0;
                for(const Fracture& ff: datagliare){

                    //bool cut = false;
                    for(unsigned int i=0; i<ff.NumVertices; i++){
                        unsigned int vert0 = i;
                        unsigned int vert1 = i+1;
                        if(i==f.NumVertices-1){
                            vert1 = i-(f.NumVertices-1);
                        }
                        if(onSegment(trace.coordTrace[0], ff.vertices[vert0], ff.vertices[vert1]) || onSegment(trace.coordTrace[1], ff.vertices[vert0], ff.vertices[vert1])){
                            vector<Fracture> nfs = cuttingfractures(ff, trace);
                            //vecCut[k++]=nfs;
                            cutted = true;
                            //newsize+=nfs.size();
                            /*if(nfs.size()==1)
                                newfractures.resize(newfractures.size()-1);*/
                            //for(const Fracture& nf: nfs){
                            //    newfractures.push_back(nf);
                            //}
                            //nfs.clear();
                            //break;
                        //}

                    //}
                    //if(!cutted){
                    //    newfractures.push_back(ff);
                    //}

                //}
                //datagliare = newfractures;
                /*unsigned int c = 0;
                unsigned int s = newfractures.size();
                for(unsigned int i=0;  i<s; i++){
                    if(newfractures[i] == Fracture{})
                        c++;
                }
                newfractures.resize(s-c);*/
                //datagliare.clear();
                //datagliare.resize(newfractures.size());
                //for(unsigned int i=0; i<newfractures.size();i++)
                //    datagliare[i]=newfractures[i];
                //newfractures.clear();
                //datagliare.resize(newfractures.size());
                //datagliare = newfractures;
                /*datagliare.resize(newsize);
                unsigned int j = 0 ;
                for(const auto& cut:vecCut){
                    for(const Fracture& fra: cut){
                        datagliare[j++] = fra;
                    }
                }*/
            //}*/


            /*for(const Fracture& nf: datagliare){
                //MapNewFractures[id++] = nf;
                cout << id << ": " << endl;
                for(unsigned int j=0; j<nf.NumVertices; j++){
                    for(unsigned int k = 0; k<3; k++){
                        cout << nf.vertices[j][k] << " ";
                    }
                    cout << endl;
                }
                id++;
            }
        }
        // dalle fratture ricavati lati e vertici -> trova modo per non mettere gloi stessi due volte

        return pm;
    }*/


    PolygonalMesh newpolygon(FractureMesh& mesh){

        PolygonalMesh pm;

        unsigned int id = 0;
        //map<unsigned int, Fracture> MapNewFractures = {};

        for(const Fracture& f: mesh.MapFractures){
            cout << "taglio della frattura " << f.id << endl;
            //Fracture f = coppia.second;

            /*map<unsigned int, Vector3d> MapCell0D = {};
            map<unsigned int, list<unsigned int>> MapCell1D = {};
            map<unsigned int, list<unsigned int>> MapCell2DEdges = {};*/

            list<Fracture> datagliare;
            datagliare.push_back(f);

            for(const Trace& trace: f.listPas){
                //cout << "la figura " << coppia.first << " ha numero tracce passanti: " << f.listPas.size() << endl;
                cout << "trace id:" << trace.id << endl;

                list<Fracture> newfractures;
                for(const Fracture& ff: datagliare){
                    vector<Fracture> nfs = cuttingfractures(ff, trace);
                    for(const Fracture& nf: nfs){
                        newfractures.push_back(nf);

                    }
                }
                datagliare.clear();
                datagliare = newfractures;
            }

            /*for(const Fracture& nf: datagliare){
                //MapNewFractures[id++] = nf;
                cout << id << ": " << endl;
                cout << nf.NumVertices<<endl;
                for(unsigned int j=0; j<nf.NumVertices; j++){
                    for(unsigned int k = 0; k<3; k++){
                        cout << nf.vertices[j][k] << " ";
                    }
                    cout << endl;
                }
                id++;
            }*/



            for(const Trace& trace: f.listNonpas){
                //cout << "la figura " << coppia.first << " ha numero tracce non passanti: " << f.listNonpas.size() << endl;
                cout << "trace (non passante) id:" << trace.id << endl;
                list<Fracture> newfractures = {};
                for(const Fracture& ff: datagliare){

                    bool cut = false;
                    if(isPointInsideFracture(trace.coordTrace[0], ff) | isPointInsideFracture(trace.coordTrace[1], ff)){
                        vector<Fracture> nfs = cuttingfractures(ff, trace);
                        cut = true;
                        for(const Fracture& nf: nfs){
                            newfractures.push_back(nf);
                        }
                    }
                    /*for(unsigned int i=0; i<ff.NumVertices; i++){
                        unsigned int &vert0 = i;
                        unsigned int vert1 = i+1;
                        if(i==f.NumVertices-1){
                            vert1 = i-(f.NumVertices-1);
                        }
                        /*for(unsigned int co =0; co<3; co++){
                            if(trace.coordTrace[0][co]>= min(ff.v)
                        }*/


                    //}
                    if(!cut){
                        //cout << "questa traccia non la taglia" << endl;
                        newfractures.push_back(ff);
                    }
                }
                datagliare = newfractures;
            }


            for(const Fracture& nf: datagliare){
                //MapNewFractures[id++] = nf;
                cout << id << ": " << endl;
                //cout << nf.NumVertices<<endl;
                for(unsigned int j=0; j<nf.NumVertices; j++){
                    for(unsigned int k = 0; k<3; k++){
                        if(abs(nf.vertices[j][k])<tol)
                            cout << scientific << setprecision(16) << 0.0 << " ";
                        else
                            cout << scientific << setprecision(16) << nf.vertices[j][k] << " ";
                    }
                    cout << endl;
                }
                id ++;
            }
            datagliare.clear();
        }
        // dalle fratture ricavati lati e vertici -> trova modo per non mettere gloi stessi due volte

        return pm;
    }


}
