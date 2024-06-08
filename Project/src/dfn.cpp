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
using namespace GeometryHelpers;

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



    vector<Fracture> cuttingfractures(Fracture& f, const Trace& t, PolygonalMesh& polyMesh){
        vector<Fracture> newfractures(2);

        vector<unsigned int> coordt = intersezionipoligonorettaLATI(t.retta, t.coordTrace[1], f);
        vector<Vector3d> coordtVere;

        if(f.tips.count(t.id)>0 && !f.tips.at(t.id))
            coordtVere = t.coordTrace;
        else
            coordtVere = intersezionipoligonoretta(t.retta, t.coordTrace[1], f);

        if(coordt.empty() || coordt.size()==1){  //la retta non interseca il poligono, aggiungo solo la fracture gia esistente
            newfractures.resize(1);
            newfractures[0] = f;
            return newfractures;
        }

        list<unsigned int> vertOldFrac = f.idvertici;
        list<unsigned int> latiOldFrac = f.idlati;

        unsigned int numNewVertici = f.NumVertices+2;
        vector<bool> angolo(2, false);
        for(const auto& vertice: f.vertices){
            for(unsigned int a = 0; a<2; a++){
                if(vertice==coordtVere[a]){
                    numNewVertici--;
                    angolo[a] = true;
                    coordt[a]--;
                }
            }
        }

        vector<Vector3d> newVertVec(numNewVertici);
        unsigned int j=0;
        unsigned int index = 0;
        //unsigned int uno = 1;
        auto it = vertOldFrac.begin();
        auto itLati = latiOldFrac.begin();
        //auto itlativecchi = f.idlati.begin();
        for(int i=0; i<numNewVertici; i++){
            if(i==coordt[index]){  // mi trovo nel punto in cui devo inserire il vertice nuovo
                if(angolo[index])  // se il vertice nuovo corrisponde ad uno che esisteva già va solo avanti con le iteraizoni
                    j++;
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
            Vector3d planeF = calcoloPiano(f);
            //double dF = -(planeF[0] * f.vertices[0][0]) - (planeF[1] * f.vertices[0][1]) - (planeF[2] * f.vertices[0][2]);
            //f.plane.reserve(4);

            //ciclo sui poigoni successivi
            for(unsigned int i=id+1; i<mesh.NumFractures; i++){
                bool inter = false;
                vector<Vector3d> trace = {};
                trace.resize(2);
                Vector3d t = {};
                double distanzainF = 0;
                double distanzainFC = 0;
                Fracture &fConf = mesh.MapFractures[i];
                f.isOnEdge = false;
                fConf.isOnEdge =false;
                Vector3d planeFConf = calcoloPiano(fConf);
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
                                if(sameLine(retta1, fConf, interFC)){
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

                else if((dist(f.barycentre, fConf.barycentre) - (maxDist(f) + maxDist(fConf))) > tol ){
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
        vector<PolygonalMesh> newfigures(mesh.NumFractures);
        //map<unsigned int, Fracture> MapNewFractures = {};

        unsigned int id = 0;
        for(Fracture& f: mesh.MapFractures){

            PolygonalMesh pm;


            cout << "taglio della frattura " << f.id << endl;
            //Fracture f = coppia.second;

            //unsigned int id0d = 0;
            //unsigned int id1d = 0;
            unsigned int id2d = 0;
            /*pm.Cell0DId.push_back(pm.numCell0D);
            pm.MapCell0D[pm.numCell0D] = f.vertices[0];
            f.idvertici.push_back(pm.numCell0D++);*/

            for(unsigned int vert = 0; vert < f.NumVertices; vert ++){
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

            list<Fracture> datagliare = {};
            datagliare.push_back(f);

            for(const Trace& trace: f.listPas){
                //cout << "la figura " << coppia.first << " ha numero tracce passanti: " << f.listPas.size() << endl;
                cout << "trace id:" << trace.id << endl;

                list<Fracture> newfractures = {};
                for(Fracture& ff: datagliare){
                    if(!f.onEdge.at(trace.id)){
                        vector<Fracture> nfs = cuttingfractures(ff, trace, pm);
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

            }



            for(const Trace& trace: f.listNonpas){
                vector<Vector3d> copiacordiTrace = trace.coordTrace;
                sort(copiacordiTrace.begin(), copiacordiTrace.end(), compareFirstElement);
                cout << "trace (non passante) id:" << trace.id << endl;
                list<Fracture> newfractures = {};
                for(Fracture& ff: datagliare){

                    //cout << "nuova ff da tyagliare" << endl;

                    bool estremiTracciaInside = false;
                    vector<Vector3d> interRetta = intersezionipoligonoretta(trace.retta, copiacordiTrace[0], ff);

                    if(interRetta.size()==2){
                        sort(interRetta.begin(), interRetta.end(), compareFirstElement);
                        vector<Vector3d> estremi(2);
                        intersezioniSuRetta(estremiTracciaInside, estremi, interRetta, copiacordiTrace);
                    }

                    if(estremiTracciaInside && !f.onEdge.at(trace.id)){
                        vector<Fracture> nfs = cuttingfractures(ff, trace, pm);
                            //cut = true;
                        for(const Fracture& nf: nfs){
                            newfractures.push_back(nf);
                        }
                    }
                    else{

                        newfractures.push_back(ff);
                    }
                }

                datagliare.clear();
                datagliare = newfractures;
            }

            for(const Fracture& nf: datagliare){
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
            }
            datagliare.clear();
            newfigures[id++] = pm;
        }

        return newfigures;
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



    void printingPolygonMesh(const vector<PolygonalMesh>& pp, const string& file){

        string filepath = "newPolygons_" + file.substr(4,file.length()-4);
        ofstream outfile(filepath);
        if (outfile.is_open()) {
            for(const PolygonalMesh& pm: pp){
                //metti tutti conrrolli per liste vuote

                outfile << "Cell0Ds: " << endl << "Id; X; Y; Z;" << endl;
                for(const auto& d: pm.Cell0DId){
                    outfile << d << "; " << pm.MapCell0D.at(d)[0] << "; " << pm.MapCell0D.at(d)[1] << "; " << pm.MapCell0D.at(d)[2] << endl;
                }
                outfile << "Cell0Ds: " << endl << "Id; IdVert1; IdVert2" << endl;
                for(const auto& e: pm.Cell1DId){
                    outfile << e << "; " << pm.MapCell0D.at(e)[0] << "; " << pm.MapCell0D.at(e)[1] << endl;
                }

                for(const auto &f: pm.Cell2DId){
                    outfile << "# FractureId; NumCell0d" << endl << f << "; " << pm.MapCell2DVertices.at(f).size() << endl << "# Cell0dId; X; Y; Z;" << endl;
                    for(const auto &t: pm.MapCell2DVertices.at(f)){
                        //auto it = pm.MapCell2DVertices.at(f).begin();
                        outfile << t << "; " ;
                        for(unsigned int k=0; k<3; k++)
                            outfile << scientific << setprecision(16)<< pm.MapCell0D.at(t)[k] << "; ";
                        outfile << endl;
                    }

                    outfile << "# Cell1ds; IdCell0ds; IdCell0ds;" << endl;
                    for(const auto &t: pm.MapCell2DEdges.at(f)){
                        //auto it = pm.MapCell2DEdges.at(f).begin();
                        outfile <<t << "; " ;
                        for(unsigned int k=0; k<2; k++)
                            outfile << scientific << setprecision(16)<< pm.MapCell1D.at(t)[k] << "; ";
                        outfile << endl;
                    }
                }
                outfile << endl;
            }
            outfile.close();
        } else {
            cout << "Impossibile aprire il file." << endl;
        }
    }






}
