#include "DFN.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>

using namespace Eigen;

namespace FractureLibrary {

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

// importo i dati dai file ( va bene per tutti i file questa funzione )
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
        mesh.FractureId.push_back(id);
        f.NumVertices = numvertices;

        cout << "la fracture  " << mesh.FractureId[i] << " ha " << numvertices << " vertici" << endl;
        getline(file,line); // leggo riga vertici senza stamparla (vediamo se farlo tutto insieme le eliminazioni)


        Vector3d bary;
        vector<Vector3d> vec;
        vec.resize(numvertices);
        for(unsigned int j=0; j<3; j++){
            getline(file, line);
            replace(line.begin(),line.end(),';',' ');
            istringstream converter(line);
            bary[i] = 0.0;

            for(unsigned int t =0; t<numvertices; t++){
                double a = 0.0;
                converter >> a;
                vec[t][j] = a;
                bary[i]+=a;
                //cout << vec[t][j];
                //listLines.push_back(line); // aggiungo alla fine della lista ogni nuova riga presa dal file
            }
            bary[j] = bary[j]/numvertices;
        }

        f.vertices = vec;
        mesh.MapFractures[id]= f;


        /*cout << mesh.MapFractures.at(id).NumVertices << endl;
        cout << (mesh.MapFractures[id]).vertices[0][0] << endl;
        for(unsigned int k=0;k<f.NumVertices;k++){
            cout << vec[k][0] << "  " << vec[k][1] << "  " << vec[k][2] << endl;
        }*/
    }

    file.close();
    return true;

}

void findIntersections(const unsigned int &id, FractureMesh &mesh){
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
        Vector4d nF;
        nF[0] = (lato1F[1]*lato2F[2] - lato1F[2]*lato2F[1])/(lato1Fn*lato2Fn);
        nF[1] = (lato1F[2]*lato2F[0] - lato1F[0]*lato2F[2])/(lato1Fn*lato2Fn);
        nF[2] = (lato1F[0]*lato2F[1] - lato1F[1]*lato2F[0])/(lato1Fn*lato2Fn);
        nF[3] = -(nF[0] * f.vertices[0][0]) - (nF[1] * f.vertices[0][1]) - (nF[2] * f.vertices[0][2]);

        for(unsigned int j =0; j<3; j++){
            lato1FConf[j] = fConf.vertices[1][j] - fConf.vertices[0][j];
            lato2FConf[j] = fConf.vertices[3][j] - fConf.vertices[0][j];
        }

        //normalizzo i lati del secondo poligono
        double lato1FnConf = lato1FConf.squaredNorm();
        double lato2FnConf = lato2FConf.squaredNorm();
        vector<double> nFConf;
        nFConf.reserve(4);
        nFConf[0] = (lato1FConf[1]*lato2FConf[2] - lato1FConf[2]*lato2FConf[1])/(lato1FnConf*lato2FnConf);
        nFConf[1] = (lato1FConf[2]*lato2FConf[0] - lato1FConf[0]*lato2FConf[2])/(lato1FnConf*lato2FnConf);
        nFConf[2] = (lato1FConf[0]*lato2FConf[1] - lato1FConf[1]*lato2FConf[0])/(lato1FnConf*lato2FnConf);
        nFConf[3] = -(nFConf[0] * fConf.vertices[0][0]) - (nFConf[1] * fConf.vertices[0][1]) - (nFConf[2] * fConf.vertices[0][2]);



        if(nF[0]==nFConf[0] && nF[1]==nFConf[1] && nF[2]==nFConf[2] && nF[3]==nFConf[3]){ //condizone di complanarità
            cout << "le figure sono complanari" << endl;
        }
        else if(nF[0]==nFConf[0] && nF[1]==nFConf[1] && nF[2]==nFConf[2] && nF[3]!=nFConf[3]){ //tolgo poligono che stanno su piani paralleli
            cout << "le figure appartengono a piani paralleli" << endl;
        }

        else if(dist(f.barycentre, fConf.barycentre) > (maxDist(f) + maxDist(fConf))){
            cout << "le figure sicuramente non si intersecano" << endl;

        }
        else{ //se non sono complanari calcolo le intersezioni tra i piani
            //Vector3d t = nF.cross(nFConf);
            cout << "bau\n";
        }




    }
}
}
