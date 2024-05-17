#include "dfn.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
#include <algorithm>

using namespace Eigen;

namespace FractureLibrary {

    VectorXd PALUSolver(const MatrixXd& a,
                    const VectorXd& b)
    {
        VectorXd solutionPALU = a.fullPivLu().solve(b);
        return solutionPALU;

    }

    bool compareFirstElement(const Vector2D& a, const Vector2D& b) {
        return a[0] < b[0];
    }


    //calcolo distanza tra due punti
    double dist(Vector3d v1, Vector3d v2){
        Vector3d d;
        for(unsigned int j =0; j<3; j++){
            d[j] = abs(v1[j] - v2[j]);
        }
        double distanza = d.squaredNorm();
        return distanza;
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

    // importo i dati dai file ( va bene per tutti questa funzione!!!!
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

        cout << "numero di fratture: " << mesh.NumFractures << endl;

        char sep;
        //list<string> listLines;
        unsigned int id = 0;
        //unsigned int numvertices = 0;
        for (i=0; i < N; i++)
        {
            Fracture f;
            getline(file,line);
            getline(file, line); //leggo l'identificatore

            istringstream converter(line);
            converter >> id >> sep >> f.NumVertices;
            f.id = id;
            mesh.FractureId[i]=id;
            //f.NumVertices = numvertices;

            cout << "la fracture  " << mesh.FractureId[i] << " ha " << numvertices << " vertici" << endl;
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

                    //listLines.push_back(line); // aggiungo alla fine della lista ogni nuova riga presa dal file
                }
                bary[j] = bary[j]/numvertices;
            }

            mesh.MapFractures[id]=f;

        }

        file.close();
        return true;

    }

    void findIntersections(const unsigned int &id, FractureMesh &mesh){
        double tol = 10 * numeric_limits<double>::epsilon();
        Fracture f = mesh.MapFractures.at(id);
        for(unsigned int i=id+1; i<mesh.NumFractures; i++){
            Fracture fConf = mesh.MapFractures.at(i);
            Vector3d lato1F;
            Vector3d lato2F;
            Vector3d lato1FConf;
            Vector3d lato2FConf;
            //calcolo 2 lati
            for(unsigned int j =0; j<3; j++){
                lato1F[j] = f.vertices[1][j] - f.vertices[0][j];
                lato2F[j] = f.vertices[3][j] - f.vertices[0][j];
            }
            //normalizzo i lati
            double lato1Fn = lato1F.squaredNorm();
            double lato2Fn = lato2F.squaredNorm();
            //calcolo il vettore normale e lo normalizzo
            Vector3d planeF;
            planeF[0] = (lato1F[1]*lato2F[2] - lato1F[2]*lato2F[1])/(lato1Fn*lato2Fn);
            planeF[1] = (lato1F[2]*lato2F[0] - lato1F[0]*lato2F[2])/(lato1Fn*lato2Fn);
            planeF[2] = (lato1F[0]*lato2F[1] - lato1F[1]*lato2F[0])/(lato1Fn*lato2Fn);
            double dF = -(planeF[0] * f.vertices[0][0]) - (planeF[1] * f.vertices[0][1]) - (planeF[2] * f.vertices[0][2]);

            for(unsigned int j =0; j<3; j++){
                lato1FConf[j] = fConf.vertices[1][j] - fConf.vertices[0][j];
                lato2FConf[j] = fConf.vertices[3][j] - fConf.vertices[0][j];
            }

            //normalizzo i lati del secondo poligono
            double lato1FnConf = lato1FConf.squaredNorm();
            double lato2FnConf = lato2FConf.squaredNorm();
            Vector3d planeFConf;
            planeFConf[0] = (lato1FConf[1]*lato2FConf[2] - lato1FConf[2]*lato2FConf[1])/(lato1FnConf*lato2FnConf);
            planeFConf[1] = (lato1FConf[2]*lato2FConf[0] - lato1FConf[0]*lato2FConf[2])/(lato1FnConf*lato2FnConf);
            planeFConf[2] = (lato1FConf[0]*lato2FConf[1] - lato1FConf[1]*lato2FConf[0])/(lato1FnConf*lato2FnConf);
            double dFConf = -(planeFConf[0] * fConf.vertices[0][0]) - (planeFConf[1] * fConf.vertices[0][1]) - (planeFConf[2] * fConf.vertices[0][2]);



            if(planeF[0]==planeFConf[0] && planeF[1]==planeFConf[1] && planeF[2]==planeFConf[2] && dF == dFConf){ //condizone di complanarità
                cout << "le figure sono complanari" << endl;
                //codice su possibile tracce lungo i lati
                //continue;
            }

            else if(planeF[0]==planeFConf[0] && planeF[1]==planeFConf[1] && planeF[2]==planeFConf[2] && dF != dFConf){
                cout << "le figure giacciono su piani paralleli" << endl;
                continue;
            }

            else if((dist(f.barycentre, fConf.barycentre) - (maxDist(f) + maxDist(fConf))) > tol ){
                cout << "le figure sicuramente non si intersecano" << endl;
                continue;

            }
            else { //se non sono complanari calcolo le intersezioni tra i piani
                Vector3d t = planeF.cross(planeFConf);
                Matrix3d A;
                A.row(0) = planeF;
                A.row(1) = planeFConf;
                A.row(2) = t;
                Vector3d b;
                b[0] = dF;
                b[1] = dFConf;
                b[2] = 0.0;

                if(A.determinant() == 0){
                    cout << "determinate = 0" << endl;
                    break;
                }

                Vector3d p = PALUSolver(A, b);

                bool onEdgeF = false;
                bool onEdgeFC = false;
                //calcolo le interseioni dei lati della fracture con la retta di interseione dei piani
                vector<Vector2d> intersectionsF;
                intersectionsF.reserve(2);
                for(unsigned int l = 0; l < f.NumVertices; l++){
                    Matrix2d A1;
                    int vert0 = l;
                    int vert1;

                    if(l<f.NumVertices-1){
                        vert1 = l +1;
                    }
                    else{
                        vert1 = l-3;
                    }
                    Vector2d r1(f.vertices[vert1][0] - f.vertices[vert0][0], p[0]+t[0] - p[0]);
                    Vector2d r2(f.vertices[vert1][1] - f.vertices[vert0][1], p[1]+t[1] - p[1]);
                    A1.row(0) = r1;
                    A1.row(1) = r2;

                    if(A1.determinant()==0){ // il lato e retta di intersezione tra i piani sono paralleli oppure sovrapposto
                        if((f.vertices[vert0][0]-p[0])/t[0] == (f.vertices[vert0][1]-p[1])/t[1]){ //questo dovrebbe dirmi se stanno sulla stessa retta
                            onEdgeF = true;
                            Vector2d v1(f.vertices[vert0][0],f.vertices[vert0][1]);
                            Vector2d v2(f.vertices[vert1][0],f.vertices[vert1][1]);
                            intersectionsF.push_back(v1);
                            intersectionsF.push_back(v2);
                            break; // esco dal ciclo for sui lati perchè ho già determinato i due punti di intersezione
                                // e non è possibile che una retta inconti un poligono conesso in più di due punti
                        }
                        else{ // le due rette sono parallele quindi non hanno intersezioni
                            continue; //dovrebbe passare al prossimo lato CONTROLLARE
                        }
                    }

                    else{
                        Vector2d b1(p[0] - f.vertices[l][0], p[1] - f.vertices[l][1]);
                        Vector2d coeff = PALUSolver(A1, b1);

                        if(coeff[0]>=0 && coeff[0]<=1){ //controllo che l'intersezione sia interna al lato del poligono (combinazione convessa)
                            Vector2d v;
                            v[0] = f.vertices[l][0] + coeff[0]*(f.vertices[l+1][0] - f.vertices[l][0]);
                            v[1] = f.vertices[l][1] + coeff[1]*(f.vertices[l+1][1] - f.vertices[l][1]);
                            intersectionsF.push_back(v);
                        }

                        if(intersectionsF.size()==2){
                            break;  //esco dal ciclo perchè di nuovo non è posisbile che una retta intersechi un poligono in più di due punti
                        }           // potremmo anche mettere questo if fuori dall'else e levare il break prima , come ci piace di più
                    }
                }

                if(intersectionsF.size()==0){ // se nessuno dei lati del poligono interseca la retta allora sicuramente i due poligono in esame
                                              // non si iintersecano -> passo a confromtarlo con un altro poligono
                    continue;
                }

                //calcolo le intersezioni della retta con l'altra fracture
                vector<Vector2d> intersectionsFC;
                intersectionsFC.reserve(2);
                for(unsigned int l = 0; l < fConf.NumVertices; l++){
                    Matrix2d A1;
                    int vert0 = l;
                    int vert1;

                    if(l<fConf.NumVertices-1){
                        vert1 = l+1;
                    }
                    else{
                        vert1 = l-3;
                    }
                    Vector2d r1(fConf.vertices[vert1][0] - fConf.vertices[vert0][0], p[0]+t[0] - p[0);
                    Vector2d r2(fConf.vertices[vert1][1] - fConf.vertices[vert0][1], p[1]+t[1] - p[1]);
                    A1.row(0) = r1;
                    A1.row(1) = r2;

                    if(A1.determinant()==0){
                        if((fConf.vertices[vert0][0]-p[0])/t[0] == (fConf.vertices[vert0][1]-p[1])/t[1]){
                            onEdgeFC = true;
                            Vector2d v1(fConf.vertices[vert0][0],fConf.vertices[vert0][1]);
                            Vector2d v2(fConf.vertices[vert1][0],fConf.vertices[vert1][1]);
                            intersectionsFConf.push_back(v1);
                            intersectionsFConf.push_back(v2);
                            break;
                        }
                        else{ //sono parallele
                            continue;
                        }
                    }

                    else{
                        Vector2d b1(p[0] - fConf.vertices[l][0], p[1] - fConf.vertices[l][1]);
                        Vector2d coeff = PALUSolver(A1, b1);

                        if(coeff[0]>=0 && coeff[0]<=1 && coeff[1] >= 0 && coeff[1]<=1){
                            Vector2d v;
                            v[0] = fConf.vertices[l][0] + coeff[0]*(fConf.vertices[l+1][0] - fConf.vertices[l][0]);
                            v[1] = fConf.vertices[l][1] + coeff[1]*(fConf.vertices[l+1][1] - fConf.vertices[l][1]);
                            intersectionsFC.push_back(v);
                        }

                        if(intersectionsFC.size()==2){
                            break;
                        }
                    }
                }

                if(intersectionsFC.size()==0){
                    continue;
                }

                // calcolo le interseioni tra i due segmenti tovati
                vector<Vector2d> trace2d;
                trace2d.resize(2);

                Matrix2d A2;
                A2.row(0)(intersectionsF[1][0] - intersectionsF[0][0], intersectionsFC[1][0] - intersectionsFC[0][0]);
                A2.row(1)(intersectionsF[1][1] - intersectionsF[0][1], intersectionsFC[1][1] - intersectionsFC[0][1]);
                Vector2d b2(intersectionsFC[0][0] - intersectionsF[0][0], intersectionsFC[0][1] - intersectionsF[0][1]);
                if(A2.determinant() == 0){ //i due segmenti stanno sulla stessa retta
                    sort(intersectionsF.begin(), intersectionsF.end(), compareFirstElement);
                    sort(intersectionsFC.begin(), intersectionsFC.end(), compareFirstElement);
                        //Vector2d coordx1(intersectionsF[0][0], intersectionsF[1][0], c
                    //coordx1.sort();
                    //Vector2d coordx2(intersectionsFC[0][0], intersectionsF[1][0]);
                    //coordx2.sort();
                    if(intersectionsF[0][0] <= intersectionsFC[0][0]){
                        if(intersectionsF[1][0]<intersectionsFC[0][0]){
                            cout << "non si intersecano" << endl;
                            continue;
                        }
                        else if(intersectionsF[1][0] >= intersectionsFC[0][0] && intersectionsF[1][0]<= intersectionsFC[1][0]){
                            trace2d(intersectionsFC[0], intersectionsF[1]);

                        }
                        else if(coordx1[1]>=coordx2[1]){
                            trace2d(intersectionsFC[0], intersectionsFC[1]);
                        }
                    }else{
                        if(intersectionsFC[1][0]<intersectionsF[0][0]){
                            cout << "non si intersecano" << endl;
                            continue;
                        }
                        else if(intersectionsFC[1][0]>=intersectionsF[0][0] && intersectionsFC[1][0]<=intersectionsF[1][0]){
                            trace2d(intersectionsF[0], intersectionsFC[1]);
                        }
                        else if(intersectionsFC[1][0]>=intersectionsF[1][0]){
                            trace2d(intersectionsF[0], intersectionsF[1]);
                        }
                    }

                    //manca controllo sulla z

                }

                else{
                    Vector2d coeff = PALUSolver(A2, b2); //controlla b2
                    if(coeff[0]>=0 && coeff[0]<=1 && coeff[1]>=0 && coeff[1]<=1){ //controllo che il punto di intersezione sia interno ai segmenti
                    Vector2d pIn(intersectionsF[0][0] + coeff[0]*(intersectionsF[1][0] - intersectionsF[0][0]),
                                     intersectionsF[0][1] + coeff[1]*(intersectionsF[1][1] - intersectionsF[0][1]));

                            //controllo che anche la coordinata z sia comune

                    }
                }
            }
        }
    }
}
