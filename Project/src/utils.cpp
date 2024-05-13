#include "dfn.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

namespace FractureLibrary {

    /*  //da sistemare per stampare i dati del nostro file...
 * bool ImportMesh(const string& filepath,
                FractureMesh& mesh)
{

    if(!ImportFR_data(filepath + "FR3_data.txt",
                       mesh))
    {
        return false;
    }
    else
    {
        cout << "C:" << endl;
        for(auto it = mesh.Cell0DMarkers.begin(); it != mesh.Cell0DMarkers.end(); it++)
        {
            cout << "key:\t" << it -> first << "\t values:";
            for(const unsigned int id : it -> second)
                cout << "\t" << id;

            cout << endl;
        }
    }*/



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

        cout << "numero di fratture: " << mesh.NumFractures << endl;

        char sep;
        list<string> listLines;
        unsigned int id = 0;
        unsigned int numvertices = 0;
        for (i=0; i < N; i++)
        {
            getline(file,line);
            getline(file, line); //leggo l'identificatore

            istringstream converter(line);
            converter >> id >> sep >> numvertices;
            mesh.FractureId.push_back(id);
            //mesh.NumVertices = numvertices;

            cout << "la fracture  " << mesh.FractureId[i] << " ha " << numvertices << " vertici" << endl;
            getline(file,line); // leggo riga vertici senza stamparla (vediamo se farlo tutto insieme le eliminazioni)

            vector<Vector3d> vec;
            vec.resize(numvertices);
            for(unsigned int j=0; j<3; j++){
                getline(file, line);
                replace(line.begin(),line.end(),';',' ');
                istringstream converter(line);

                for(unsigned int i =0; i<numvertices; i++){
                    double a = 0.0;
                    converter >> a;
                    vec[i][j] = a;

                    listLines.push_back(line); // aggiungo alla fine della lista ogni nuova riga presa dal file
                }
            }

            mesh.MapFractures.insert({id,{vec}});

            //calcolo le intersezioni tre le figure
            /*for(unsigned int i = 0; i<mesh.NumFractures -1; i++){
                for(unsigned int j = 1; j < mesh.NumFractures; j++){
                    //confronto ogi figura con tutte quelle che vengono dopo di lei
                    Vector3d p;
                    for(unsigned int k = 0; k < mesh.MapFractures.at(i).size(); k++){
                        p[1]= mesh.MapFractures.at(i)[k][1];
                        p[2]= mesh.MapFractures.at(i)[k][2];
                        p[3]= mesh.MapFractures.at(i)[k][3];
                    }

                }
            }*/

            //for(unsigned int i = 0; i < numvertices; i++)
            //  converter >> x[i];

            /*getline(file, line);
            replace(line.begin(),line.end(),';',' ');

            listLines.push_back(line); // aggiungo alla fine della lista ogni nuova riga presa dal file
             // lista coordinate y
            //for(unsigned int i = 0; i < numvertices; i++)
            //  converter >> y[i];

            getline(file, line);
            replace(line.begin(),line.end(),';',' ');
            listLines.push_back(line); // aggiungo alla fine della lista ogni nuova riga presa dal file
             // lista coordinate z
            //for(unsigned int i = 0; i < numvertices; i++)
            //  converter >> z[i];
            vector<double> v(3);

            for (unsigned int i = 0; i < numvertices; i++)
                converter >> x[i] >> y[i] >> z[i];
            v[0] = x[i];
            v[1] = y[i];
            v[2] = z[i];
            mesh.CoordinatesFractures(v);


            // Skip Comment Line
            getline(file,line);*/

        }



        // qui comincio a dividere ogni riga in id, marker e coordinate (x,y)


        file.close();
        return true;
    }
}
