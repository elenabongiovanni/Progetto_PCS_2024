#include "DFN.hpp"
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
    string Id;
    list<string> listLines;
    unsigned int id = 0;
    unsigned int numvertices = 0;
    for (i=0; i< N; i++)
    {
        getline(file,Id, ';'); //leggo l'identificatore
        getline(file,line);
        istringstream converter(line);
        converter >> id >> numvertices;
        mesh.FractureId.push_back(id);
        mesh.NumVertices = numvertices;
        getline(file,line); // leggo riga vertici senza stamparla

        getline(file, line);
        replace(line.begin(),line.end(),';',' ');
        listLines.push_back(line); // aggiungo alla fine della lista ogni nuova riga presa dal file
        vector<double> x(numvertices); // lista coordinate x
        //for(unsigned int i = 0; i < numvertices; i++)
          //  converter >> x[i];

        getline(file, line);
        replace(line.begin(),line.end(),';',' ');
        listLines.push_back(line); // aggiungo alla fine della lista ogni nuova riga presa dal file
        vector<double> y(numvertices); // lista coordinate y
        //for(unsigned int i = 0; i < numvertices; i++)
          //  converter >> y[i];

        getline(file, line);
        replace(line.begin(),line.end(),';',' ');
        listLines.push_back(line); // aggiungo alla fine della lista ogni nuova riga presa dal file
        vector<double> z(numvertices); // lista coordinate z
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
        getline(file,line);

    }



    // qui comincio a dividere ogni riga in id, marker e coordinate (x,y)


    file.close();
    return true;
}
}