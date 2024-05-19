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

    bool compareFirstElement(const Vector3d& a, const Vector3d& b) {
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

            f.barycentre = bary;
            mesh.MapFractures[id]=f;

        }

        file.close();
        return true;

    }

    void findIntersections(const unsigned int &id, FractureMesh &mesh){
        double tol = 10 * numeric_limits<double>::epsilon();
        mesh.vecTrace.reserve(mesh.NumFractures*(mesh.NumFractures-1));
        Fracture f = mesh.MapFractures.at(id);
        Vector3d lato1F;
        Vector3d lato2F;
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

        //ciclo sui poigoni successivi
        for(unsigned int i=id+1; i<mesh.NumFractures; i++){
            bool inter = false;
            vector<Vector3d> trace;
            trace.resize(2);

            Fracture fConf = mesh.MapFractures.at(i);

            Vector3d lato1FConf;
            Vector3d lato2FConf;

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

            //cout <<"Il piano "<< id << " ha equazione "<< planeF[0] <<"x + " <<planeF[1] <<"y + " <<planeF[2] <<"z + " << dF << " = 0 "<<endl;
            //cout <<"Il piano "<< i << " ha equazione "<< planeFConf[0] <<"x + " <<planeFConf[1] <<"y + " <<planeFConf[2] <<"z + " << dFConf << "=0 "<<endl;

            if(planeF[0]==planeFConf[0] && planeF[1]==planeFConf[1] && planeF[2]==planeFConf[2] && dF == dFConf){ //condizone di complanarità
                cout << "le figure "<< id <<" e " << i <<" sono complanari" << endl;
                //codice su possibile tracce lungo i lati


                //continue;
            }

            else if(planeF[0]==planeFConf[0] && planeF[1]==planeFConf[1] && planeF[2]==planeFConf[2] && dF != dFConf){
                cout << "le figure "<< id <<" e " << i <<" giacciono su piani paralleli" << endl;
                continue;
            }

            else if((dist(f.barycentre, fConf.barycentre) - (maxDist(f) + maxDist(fConf))) > tol ){
                cout << "le figure "<< id <<" e " << i <<" sicuramente non si intersecano" << endl;
                continue;

            }
            else { //se non sono complanari calcolo le intersezioni tra i piani
                Vector3d t = planeF.cross(planeFConf);
                //cout << "retta di intersezione: " << t[0] << " " <<t[1] << " " << t[2] << endl;
                Matrix3d A;
                A.row(0) = planeF;
                A.row(1) = planeFConf;
                A.row(2) = t;
                Vector3d b;
                b[0] = -dF;
                b[1] = -dFConf;
                b[2] = 0.0;

                if(A.determinant() == 0){
                    //cout << "Le figure "<< id <<" e " << i <<" il deteminante è nullo ora non mi ricordo il perchè" << endl;
                    break; //questo non è di nuovo il caso di complanarità???
                }

                Vector3d p = PALUSolver(A, b);
                //cout << "punto p di intersezione: " << p[0] << " " <<p[1] << " " << p[2] << endl;

                //bool onEdgeF = false;
                //bool onEdgeFC = false;

                //calcolo le interseioni dei lati della fracture con la retta di interseione dei piani
                vector<Vector3d> intersectionsF;
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
                    Vector2d r1(f.vertices[vert1][0] - f.vertices[vert0][0], t[0]);
                    Vector2d r2(f.vertices[vert1][1] - f.vertices[vert0][1], t[1]);
                    A1.row(0) = r1;
                    A1.row(1) = r2;

                    if(A1.determinant()==0){ // il lato e retta di intersezione tra i piani sono paralleli oppure sovrapposto
                        if((f.vertices[vert0][0]-p[0])/t[0] == (f.vertices[vert0][1]-p[1])/t[1]){ //questo dovrebbe dirmi se stanno sulla stessa retta
                            //onEdgeF = true;
                            Vector3d v1(f.vertices[vert0][0],f.vertices[vert0][1], f.vertices[vert0][2]);
                            Vector3d v2(f.vertices[vert1][0],f.vertices[vert1][1], f.vertices[vert1][2]);
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
                        Vector2d b1(p[0] - f.vertices[vert0][0], p[1] - f.vertices[vert0][1]);
                        Vector2d coeff = PALUSolver(A1, b1);

                        if(coeff[0]>=0 && coeff[0]<=1){ //controllo che l'intersezione sia interna al lato del poligono (combinazione convessa)
                            //cout<<"entra qui"<<endl;
                            Vector3d v;
                            v[0] = f.vertices[vert0][0] + coeff[0]*(f.vertices[vert1][0] - f.vertices[vert0][0]);
                            v[1] = f.vertices[vert0][1] + coeff[0]*(f.vertices[vert1][1] - f.vertices[vert0][1]);
                            v[2] = (-dF - planeF[0]*v[0] - planeF[1]*v[1])/planeF[2];
                            intersectionsF.push_back(v);
                        }

                        if(intersectionsF.size()==2){
                            //cout << "il primo poligono è tagliato" << endl;
                            break;  //esco dal ciclo perchè di nuovo non è posisbile che una retta intersechi un poligono in più di due punti
                        }           // potremmo anche mettere questo if fuori dall'else e levare il break prima , come ci piace di più
                    }
                }

                if(intersectionsF.size()==0){ // se nessuno dei lati del poligono interseca la retta allora sicuramente i due poligono in esame
                                              // non si iintersecano -> passo a confromtarlo con un altro poligono
                    cout << "Le figure " << id <<" e " << i <<" non si intersecano perche' la prima non e' attraversata dalla retta di intersezione tra i piani" << endl;
                    continue;
                }


                //calcolo le intersezioni della retta con l'altra fracture
                vector<Vector3d> intersectionsFC;
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
                    Vector2d r1(fConf.vertices[vert1][0] - fConf.vertices[vert0][0], t[0]);
                    Vector2d r2(fConf.vertices[vert1][1] - fConf.vertices[vert0][1], t[1]);
                    //cout << r1[0] << " " << r2[0] << endl;
                    A1.row(0) = r1;
                    A1.row(1) = r2;

                    if(A1.determinant()==0){
                        //cout << "det = 0";
                        if((fConf.vertices[vert0][0]-p[0])/t[0] == (fConf.vertices[vert0][1]-p[1])/t[1]){
                            //onEdgeFC = true;
                            Vector3d v1(fConf.vertices[vert0][0],fConf.vertices[vert0][1], fConf.vertices[vert0][2]);
                            Vector3d v2(fConf.vertices[vert1][0],fConf.vertices[vert1][1], fConf.vertices[vert1][2]);
                            intersectionsFC.push_back(v1);
                            intersectionsFC.push_back(v2);
                            break;
                        }
                        else{ //sono parallele
                            continue;
                        }
                    }

                    else{
                        Vector2d b1(p[0] - fConf.vertices[vert0][0], p[1] - fConf.vertices[vert0][1]);
                        Vector2d coeff = PALUSolver(A1, b1);

                        if(coeff[0]>=0 && coeff[0]<=1){
                            Vector3d v;
                            v[0] = fConf.vertices[vert0][0] + coeff[0]*(fConf.vertices[vert1][0] - fConf.vertices[vert0][0]);
                            v[1] = fConf.vertices[vert0][1] + coeff[0]*(fConf.vertices[vert1][1] - fConf.vertices[vert0][1]);
                            v[2] = (- dFConf - planeFConf[0]*v[0] - planeFConf[1]*v[1])/planeFConf[2];
                            intersectionsFC.push_back(v);
                        }

                        if(intersectionsFC.size()==2){
                            //cout << "il secondo poligono è tagliato" << endl;
                            break;
                        }
                    }
                }

                if(intersectionsFC.size()==0){
                    cout << "Le figure " << id <<" e " << i << " non si intersecano perche' la seconda non e' attraversata dalla retta di intersezione tra i piani" << endl;
                    continue;
                }

                // calcolo le interseioni tra i due segmenti tovati

               // Matrix2d A2;
               // A2.row(0)(intersectionsF[1][0] - intersectionsF[0][0], intersectionsFC[1][0] - intersectionsFC[0][0]);
                //A2.row(1)(intersectionsF[1][1] - intersectionsF[0][1], intersectionsFC[1][1] - intersectionsFC[0][1]);
                //Vector2d b2(intersectionsFC[0][0] - intersectionsF[0][0], intersectionsFC[0][1] - intersectionsF[0][1]);
                //if(A2.determinant() == 0){ //i due segmenti stanno sulla stessa retta
                sort(intersectionsF.begin(), intersectionsF.end(), compareFirstElement);
                sort(intersectionsFC.begin(), intersectionsFC.end(), compareFirstElement);
                    //cout <<"intersezioni x primo poligono " << intersectionsF[0][0] << " " << intersectionsF[1][0] << endl;;
                    //cout <<"intersezionix secondo poligono " << intersectionsFC[0][0] << " " << intersectionsFC[1][0] << endl;

                if(intersectionsF[0][0] <= intersectionsFC[0][0]){
                    if(intersectionsF[1][0]<intersectionsFC[0][0]){
                        cout <<"Le figure "<< id <<" e " << i << " non si intersecano" << endl;
                        continue;
                    }
                    else if(intersectionsF[1][0] >= intersectionsFC[0][0] && intersectionsF[1][0]<= intersectionsFC[1][0]){
                        double z1F = (-dF - planeF[0]*intersectionsFC[0][0] - planeF[1]*intersectionsFC[0][1])/planeF[2];
                        double z2FC = (-dFConf - planeFConf[0]*intersectionsF[1][0] - planeFConf[1]*intersectionsF[1][1])/planeFConf[2];
                        if((intersectionsFC[0][2]-z1F)<tol && (intersectionsF[1][2]-z2FC)<tol){
                            trace.push_back(intersectionsFC[0]);
                            trace.push_back(intersectionsF[1]);
                            cout << "LE DUE FIGURE " << id <<" e " << i << " SI INTERSECANO" << endl;
                            inter = true;


                            /*if((intersectionsF[1][0]-intersectionsFC[1][0])<tol){
                                cout << "la traccia è passante per la figura " << i << endl;
                            }
                            if((intersectionsFC[0][0] - intersectionsF[0][0])<tol){
                                cout << "la traccia è passante per la figura " << id << endl;
                            }*/

                        }
                        else{
                            cout << id <<" e " << i << " in realtà non si intersecano sulla z" << endl;
                            continue;
                        }

                    }
                    else if(intersectionsF[1][0]>=intersectionsFC[1][0]){
                        double z1F = (-dF - planeF[0]*intersectionsFC[0][0] - planeF[1]*intersectionsFC[0][1])/planeF[2];
                        double z2F = (-dF - planeF[0]*intersectionsFC[1][0] - planeF[1]*intersectionsFC[1][1])/planeF[2];
                        if((intersectionsFC[0][2]-z1F)<tol && (intersectionsFC[1][2]-z2F)<tol){
                            trace.push_back(intersectionsFC[0]);
                            trace.push_back(intersectionsFC[1]);
                            cout << "LE DUE FIGURE " << id <<" e " << i <<" SI INTERSECANO" << endl;
                            inter = true;
                            //cout << "la traccia e' passante per la figura " << i << endl;
                        }
                        else{
                            cout << id <<" e " << i  << " in realtà non si intersecano sulla z" << endl;
                            continue;
                        }
                    }

                }else{
                    if(intersectionsFC[1][0]<intersectionsF[0][0]){
                        cout <<"Le figure "<< id <<" e " << i << " non si intersecano" << endl;
                        continue;
                    }
                    else if(intersectionsFC[1][0]>=intersectionsF[0][0] && intersectionsFC[1][0]<=intersectionsF[1][0]){
                        double z1FC = (-dFConf - planeFConf[0]*intersectionsF[0][0] - planeFConf[1]*intersectionsF[0][1])/planeFConf[2];
                        double z2F = (-dF - planeF[0]*intersectionsFC[1][0] - planeF[1]*intersectionsFC[1][1])/planeF[2];
                        if((intersectionsF[0][2]-z1FC)<tol && (intersectionsFC[1][2]-z2F)<tol){
                            trace.push_back(intersectionsF[0]);
                            trace.push_back(intersectionsFC[1]);
                            cout << "LE DUE FIGURE " << id <<" e " << i << " SI INTERSECANO" << endl;
                            inter = true;

                           // if()
                        }
                        else{
                            cout << id <<" e " << i  << " in realtà non si intersecano sulla z" << endl;
                            continue;
                        }

                    }
                    else if(intersectionsFC[1][0]>=intersectionsF[1][0]){
                        double z1FC = (-dFConf - planeFConf[0]*intersectionsF[0][0] - planeFConf[1]*intersectionsF[0][1])/planeFConf[2];
                        double z2FC = (-dFConf - planeFConf[0]*intersectionsF[1][0] - planeFConf[1]*intersectionsF[1][1])/planeFConf[2];
                        if((intersectionsF[0][2]-z1FC)<tol && (intersectionsF[1][2]-z2FC)<tol){
                            trace.push_back(intersectionsF[0]);
                            trace.push_back(intersectionsF[1]);
                            cout << "LE DUE FIGURE " << id <<" e " << i << " SI INTERSECANO" << endl;
                            inter = true;
                            //cout << "la traccia e' passante per la figura " << i << endl;
                        }
                        else{
                            cout << id <<" e " << i  << " in realtà non si intersecano sulla z" << endl;
                            continue;
                        }
                    }
                }

            }

            if(inter){
                Trace newTrace;
                newTrace.coordTrace.resize(2);
                newTrace.coordTrace = trace;
                mesh.vecTrace.push_back(newTrace);

                vector<bool> interLatiF(2,false);
                //interLatiF.resize(2);
                //interLatiF(false,false);

                vector<bool> interLatiFC(2,false);
                //interLatiFC.resize(2);
                //interLatiFC(false,false);


                //ciclo su tutti i lati di entrambi i poligoni per capire se gli estrami della traccia appartengono ai lati, ovvero se la traccia è passante
                for(unsigned int k = 0; k<2; k++){
                    Vector3d punto = newTrace.coordTrace[k];
                    for(unsigned int j = 0; j<4; j++){
                        Vector3d rettaF;
                        Vector3d rettaFC;
                        int vert0 = j;
                        int vert1;
                        if(j<3){
                            vert1 = j+1;
                        }
                        else{
                            vert1 = j-3;
                        }
                        rettaF[0]=f.vertices[vert1][0]-f.vertices[vert0][0];
                        rettaF[1]=f.vertices[vert1][1]-f.vertices[vert0][1];
                        rettaF[2]=f.vertices[vert1][2]-f.vertices[vert0][2];
                        rettaFC[0]=fConf.vertices[vert1][0]-fConf.vertices[vert0][0];
                        rettaFC[1]=fConf.vertices[vert1][1]-fConf.vertices[vert0][1];
                        rettaFC[2]=fConf.vertices[vert1][2]-fConf.vertices[vert0][2];


                        if(punto[0]*rettaF[0] + punto[1]*rettaF[1] + punto[2]*rettaF[2]<tol){
                            if(punto[0]>=min(f.vertices[vert1][0],f.vertices[vert0][0]) - tol && punto[1]>=min(f.vertices[vert1][1],f.vertices[vert0][1])-tol && punto[2]>=min(f.vertices[vert1][2],f.vertices[vert0][2])-tol &&
                                punto[0]<=max(f.vertices[vert1][0],f.vertices[vert0][0])+tol && punto[1]<=max(f.vertices[vert1][1],f.vertices[vert0][1])+tol && punto[2]<=max(f.vertices[vert1][2],f.vertices[vert0][2])+tol){
                                interLatiF[k] = true;
                            }
                        }
                        if(punto[0]*rettaFC[0] + punto[1]*rettaFC[1] + punto[2]*rettaFC[2]<tol){
                            if(punto[0]>=min(fConf.vertices[vert0][0],fConf.vertices[vert1][0]) - tol && punto[1]>=min(fConf.vertices[vert0][1],fConf.vertices[vert1][1])-tol && punto[2]>=min(fConf.vertices[vert0][2],fConf.vertices[vert1][2])-tol &&
                                punto[0]<=max(fConf.vertices[vert1][0],fConf.vertices[vert0][0])+tol && punto[1]<=max(fConf.vertices[vert1][1],fConf.vertices[vert0][1])+tol && punto[2]<=max(fConf.vertices[vert1][2],fConf.vertices[vert0][2])+tol){
                                interLatiFC[k] = true;
                            }
                        }
                    }
                }

                if(interLatiF[0] && interLatiF[1]){
                    newTrace.fracturesTrace.insert({id, true});
                    cout << "figura " << id << " e' passante"<< endl;
                }
                else{
                    newTrace.fracturesTrace.insert({id, false});
                }

                if(interLatiFC[0] && interLatiFC[1]){
                    newTrace.fracturesTrace.insert({i, true});
                    cout << "figura " << i << " e' passante"<<endl;
                }
                else{
                    newTrace.fracturesTrace.insert({i, false});
                }

            }
        }

    }
}
