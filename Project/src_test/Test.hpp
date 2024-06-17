#ifndef TEST_HPP
#define TEST_HPP

#include <gtest/gtest.h>
#include "DFN.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include <fstream>
#include "UCDUtilities.hpp"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

namespace FractureLibrary {


// Test per la funzione double dist(Vector3d v1, Vector3d v2);
TEST(FractureTEST, TestFunctionDistance){
    Vector3d v1 = Vector3d::Zero();
    Vector3d v2 = Vector3d::Zero();

    v1 << 3.0, sqrt(2), 2.0;
    v2 << 4.0, 4*sqrt(2), -7.0;

    double distanza = dist(v1,v2);
    double distanzaaspettata = 10.0;
    EXPECT_NEAR(distanza, distanzaaspettata, 1e-8);
}

//Test per la funzione bool orderLen(const Trace &a, const Trace &b)
TEST(FractureTEST, TestFunctionOrderLen){
    Trace t1;
    Trace t2;
    t1.len = 3.0;
    t2.len = 1.0;
    bool lunghezza = orderLen(t1,t2);
    EXPECT_EQ(lunghezza, true);
}


//Test per la funzione bool compareFirstElement(const Vector3d& a, const Vector3d& b)
TEST(FractureTEST, TestFunctioncompareFirstElement){
    Vector3d a = Vector3d::Zero();
    Vector3d b = Vector3d::Zero();

    a << 0.0, 2.0, 4.0;
    b << 3.0, 1.0, 0.0;

    bool compare = compareFirstElement(a, b); // mi dice se sono già in ordine

    EXPECT_EQ(compare, true);
}

//Test per la funzione double maxDist(Fracture f)
TEST(FractureTEST, TestFunctionmaxDist){
    Fracture f;
    f.NumVertices = 4;
    f.barycentre = {0.5, 1.0, 2.0};
    f.vertices = {{0.5, 1.0, 3.0},{0.5, 5.0, 0.0},{1.0, -2.0, 3.0},{0.0, 0.0, 2.0}};

    double max = f.maxDist();
    double maxaspettata = 2*sqrt(5);
    EXPECT_NEAR(max, maxaspettata,1e-8);
}

//Test per la funzione bool onSegment(const Vector3d& p, const Vector3d& a, const Vector3d& b)
TEST(FractureTEST, TestOnSegment){
    Vector3d p = Vector3d::Zero();
    Vector3d a = Vector3d::Zero();
    Vector3d b = Vector3d::Zero();

    p << 0.5, 1.0, 0.5;
    a << 3.0, 6.0, 3.0;
    b << 1.0, 2.0, 1.0;

    bool OnSegment = onSegment(p,a,b);

    EXPECT_TRUE(OnSegment);
}


//Test per la funzione Vector3d calcoloPiano(const Fracture& f)
TEST(FractureTEST, Testcalcolopiano){

    Fracture f;
    f.vertices = {{1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{0.0, 1.0, 0.0},{0.0, 0.0, 0.0}}; //num vertici uguale a 4

    Vector3d piano = f.calcoloPiano();
    Vector3d pianoaspettato(0.0, 0.0, 1.0);

    EXPECT_TRUE(piano.isApprox(pianoaspettato,1e-8));
}

//Test per la funzione vector<Vector3d> intersezionipoligonoretta(const Vector3d& t, const Vector3d p, const Fracture& f)
// 1° test: caso in cui la retta appartiene al piano: due punti di intersezione
TEST(FractureTEST, Testintersezionipoligonoretta2) {
    Vector3d t, p;
    p << 0.0, 1.0, 2.5;
    //t << 0.0, 0.0, -10.0;
    bool onEdge = true;
    // piano −10x+12y−14=0
    vector<Vector3d>  f = {{1.0, 2.0, 3.0}, {-1.0, 0.0, 2.0}, {-1.0, 0.0, 7.0}, {1.0, 2.0, 7.0}}; // num vertici uguale a 4
    Vector3d piano (-10,10,0);
    // piano x=0
    Vector3d pianoConf (1,0,0);
    t = piano.cross(pianoConf);
    //Vector2d beta;
    vector<Vector3d> intersezioni = intersezionipoligonoretta(t,p,f,onEdge);

    vector<Vector3d> intersezioniaspettate = {Vector3d(0.0, 1.0, 2.5),Vector3d(0.0,1.0,7.0)};

    ASSERT_EQ(intersezioni.size(), intersezioniaspettate.size());
    for (size_t i = 0; i < intersezioniaspettate.size(); ++i) {
        EXPECT_TRUE(intersezioni[i].isApprox(intersezioniaspettate[i], 1e-2));
    }
}


//Test per la funzione vector<Vector3d> intersezionipoligonoretta(const Vector3d& t, const Vector3d p, const Fracture& f)
// 2° test: caso in cui la retta non appartiene al piano: un punto di intersezione
TEST(FractureTEST, Testintersezionipoligonoretta1) {
    Vector3d t, p;
    p << 2.0, 8.0, 5.0;
    //t << -2.0, -8.0, -1.0;
    bool onEdge = false;
    Vector3d piano (4.0,-1.0,0.0);
    Vector3d pianoConf (-1.0,0.0,2.0);
    t = piano.cross(pianoConf);

    vector<Vector3d> f = {{-6.0, 5.0, 0.0}, {-6.0, 0.5, 0.0}, {0.0, 0.0, 4.0}, {0.0, 5.0, 4.0}}; // num vertici uguale a 4

    //Vector2d beta;
    vector<Vector3d> intersezioni = intersezionipoligonoretta(t, p, f,onEdge);
    vector<Vector3d> intersezioniaspettate = {Vector3d(0.0, 0.0, 4.0)};

    ASSERT_EQ(intersezioni.size(), intersezioniaspettate.size());
    for (size_t i = 0; i < intersezioniaspettate.size(); ++i) {
        EXPECT_TRUE(intersezioni[i].isApprox(intersezioniaspettate[i], 1e-2));
    }
}


//Test per la funzione VectorXd PALUSolver(const MatrixXd& a, const VectorXd& b)
TEST(FRACTURETest, PALUSolver) {
    // Definizione della matrice e del vettore
    MatrixXd A(3, 3);
    A << 2, -1, 0,
        1, 1, 1,
        0, 0, 2;

    VectorXd b(3);
    b << 0, 4, 2;

    // Calcolo della soluzione
    VectorXd solution = PALUSolver(A, b);

    // Verifica della soluzione attesa
    VectorXd expected_solution(3);
    expected_solution << 1.0, 2.0, 1.0;

    // Confronto tra la soluzione calcolata e quella attesa
    for (int i = 0; i < solution.size(); ++i) {
        EXPECT_NEAR(solution(i), expected_solution(i), 1e-8); // Utilizzo una tolleranza di 1e-6 per confrontare i numeri floating point
    }
}

//Test per la funzione bool sameLine(const Vector3d& retta, const Vector3d& p, const vector<Vector3d>& f, vector<Vector3d> coordinate);
TEST(FRACTURETest, testsameLine) {
    // retta y = 4/3x + 2/3
    Vector3d ret(3.0,4.0,0.0);
    Vector3d p(1.0,2.0,0.0);

    // dovrebbe stare sulla stessa retta con le coordinate (-0.5,0.0,0.0) e (2.5,4.0,0.0)
    vector<Vector3d> f = {{0.0,-1.0,1.0},{-0.5,0.0,0.0},{2.5,4.0,0.0},{4.5,0.0,-2.0}};
    vector<Vector3d> coordinate;

    bool sameline = sameLine(ret,p,f,coordinate);
    EXPECT_TRUE(sameline);
}

// TEST per la funzione bool intersLato(const Trace& tra, const Cell0d& c1, const Cell0d& c2, Vector3d& inters, PolygonalMesh& pm);
TEST(FRACTURETest, testintersLato) {
    PolygonalMesh pm;
    Cell0d c1;
    Cell0d c2;
    Trace t;
    // calcola le intersezioni tra la retta della traccia e le coordinate dei vertici
    Vector3d inters;
    t.retta = {1.0,1.0,0.0};
    t.coordTrace = {{0.0,0.0,1.0},{1.0,1.0,1.0}};

    c1.coordinates = {1.0,2.5,1.0};
    c2.coordinates = {4.0,2.0,1.0};

    bool interslato = intersLato(t,c1,c2,inters,pm);
    EXPECT_TRUE(intersLato);
    EXPECT_NEAR(inters[0],2.2857,1e-4);
    EXPECT_NEAR(inters[1],2.2857,1e-4);
    EXPECT_NEAR(inters[2],1.0,1e-4);


}

// TEST per la funzione bool checkIsNew(const Cell1d& c2d, const Vector3d& point, const PolygonalMesh& pm, unsigned int& id);
/*TEST(FRACTURETest, testcheckisnew) {
    PolygonalMesh pm;
    unsigned int id = 0;
    Vector3d point(1.0,2.0,2.0); //diverse!!!
    Cell1d c2d;
    unsigned int id0 = 0;
    unsigned int id1 = 1;
    Vector3d coord0 = {1.0,2.0,4.0};
    Cell0d c0d(id0,coord0);
    Vector3d coord1 = {1.0,3.0,4.0};
    Cell0d c1d(id1,coord1);
    c2d.extremes = {0,1};
    pm.MapCell0D[0] = c0d;
    pm.MapCell0D[1] = c1d;

    // verifica che le coordinate nuove inserite siano diverse da quelle già all'interno della mappa
    bool checkisnew = checkIsNew(c2d,point,pm,id);

    EXPECT_TRUE(checkisnew);
}*/

//Test per la funzione bool ImportFR_data(const string &filename, FractureMesh& mesh)
TEST(FRACTURETEST, TestPlotImportData){
    FractureMesh mesh;
    string filename = "import_test";
    ofstream outfile(filename);
    if (outfile.is_open()) {
        outfile << "# N Number of Fractures"<<endl <<"2"<<endl;
        outfile<< "# FractureId; NumVertices"<<endl<<"0;3"<<endl;
        outfile<<"#Vertices"<<endl<<  "0.0000000000000000e+00; 1.0000000000000000e+00; 1.0000000000000000e+00"
                <<endl <<  "0.0000000000000000e+00; 0.0000000000000000e+00; 1.0000000000000000e+00"
                <<endl << "0.0000000000000000e+00; 0.0000000000000000e+00; 0.0000000000000000e+00"<<endl;

        outfile<< "# FractureId; NumVertices"<<endl<<"0;3"<<endl;
        outfile<<"#Vertices"<<endl<<  "0.0000000000000000e+00; 1.0000000000000000e+00; 1.0000000000000000e+00"
                <<endl <<  "0.0000000000000000e+00; 0.0000000000000000e+00; 1.0000000000000000e+00"
                <<endl << "0.0000000000000000e+00; 0.0000000000000000e+00; 0.0000000000000000e+00"<<endl;

        outfile.close();
    }
    else {
        cout << "Impossibile aprire il file." << endl;
    }
    bool import = ImportFR_data(filename,mesh);
    EXPECT_TRUE(import);

}

//TEST per la funzione bool intersezioniSuRetta(vector<Vector3d>& trace, vector<Vector3d>& s1, vector<Vector3d>& s2)

// primo test: caso in cui si intersecano
TEST(FRACTURETEST, Testintersezionisuretta1){
    bool bole = false;
    vector<Vector3d> trace(2);

    vector<Vector3d> s1 = {{2.0,0.0,0.0},{4.5,0.0,0.0}};
    vector<Vector3d> s2 = {{4.0,0.0,0.0},{5.0,0.0,0.0}};

    // calcola le coordinate della traccia
    intersezioniSuRetta(bole,trace,s1,s2);

    EXPECT_TRUE(bole);
    EXPECT_EQ(trace[0],s2[0]);
    EXPECT_EQ(trace[1],s1[1]);

}

// secondo test: caso in cui non si intersecano
TEST(FRACTURETEST, Testintersezionisuretta2){
    bool bole = false;
    vector<Vector3d> trace(2);

    vector<Vector3d> s1 = {{2.0,0.0,0.0},{3.0,0.0,0.0}};
    vector<Vector3d> s2 = {{4.0,0.0,0.0},{5.0,0.0,0.0}};


    intersezioniSuRetta(bole,trace,s1,s2);
    EXPECT_FALSE(bole);

}


/*
//TEST per la funzione  void printingtraces(const string& file)
TEST(FRACTURETEST, Testprintingtraces){
    FractureMesh mesh;
    Trace t;
    t.fraId = {0,1};
    //map<unsigned int, bool> fracturesTrace ={};
    t.coordTrace = {{0.0,0.0,0.0},{2.0,1.0,3.0}};
    t.len = dist(t.coordTrace[0],t.coordTrace[1]);
    unsigned int id = 0;
    t.retta = {2.0,1.0,3.0};
    t.p = {0.0,0.0,0.0};
    mesh.MapTrace[0] = t;
    string file = "    test_traccia    ";
    mesh.printingtraces(file);
}
*/

//TEST per la funzione void convertFracture(const Fracture& f, unsigned int& idVert, unsigned int& idEdge, unsigned int& id2d)
TEST(FRACTURETEST, TestconvertFracture){
    Fracture f;
    f.id = 0;
    f.NumVertices = 4;
    f.vertices = {{0.0,1.0,2.0},{-4.0,2.0,2.0},{-3.0,-1.0,2.0},{0.0,-1.0,2.0}};
    Cell2d c;
    unsigned int idVert;
    unsigned int idEdge;
    unsigned int id2d;

    c.convertFracture(f,idVert, idEdge, id2d);

    // controllo che le coordinate inserite siano corrette
    for (unsigned int i = 0; i<f.NumVertices;i++){
         EXPECT_EQ(c.Cell2DVertices[i].coordinates,f.vertices[i]);
    }

    // controllo correttezza estremi dei lati
    for (unsigned int i = 0; i<f.NumVertices;i++){
         unsigned int j=(i+1)%f.NumVertices;
        EXPECT_EQ(c.Cell2DEdges[i].extremes[0],c.Cell2DVertices[i].id);
        EXPECT_EQ(c.Cell2DEdges[i].extremes[1],c.Cell2DVertices[j].id);
    }

}


//TEST per la funzione void addFirstCell2d(Cell2d& c2)
TEST(FRACTURETEST, TestaddfirstCell2d){
    PolygonalMesh pm;
    Fracture f;
    f.id = 0;
    f.NumVertices = 4;
    f.vertices = {{0.0,1.0,2.0},{-4.0,2.0,2.0},{-3.0,-1.0,2.0},{0.0,-1.0,2.0}};
    Cell2d c;
    unsigned int idVert;
    unsigned int idEdge;
    unsigned int id2d;

    c.convertFracture(f,idVert, idEdge, id2d);

    pm.addFirstCell2d(c);

    EXPECT_EQ(pm.numCell0D,4);
    //EXPECT_TRUE(pm.MapCell0D[0].old);
}


//TEST funzione void findIntersections(FractureMesh &mesh)
TEST(FRACTURETEST, Testfindintersections1){
    FractureMesh mesh;
    Fracture f1;
    Fracture f2;
    Fracture f3;
    unsigned int id0 = 0;
    unsigned int id1 = 1;
    unsigned int id2 = 2;

    f1.id = id0;
    f1.vertices = {{1.0, 2.0, 3.0}, {1.0, 2.0, 7.0}, {-1.0, 0.0, 7.0}, {-1.0, 0.0, 2.0}};
    f1.NumVertices = 4;
    f1.barycentre = {0.0,1.0,4.75};

    f2.id = id1;
    f2.vertices = {{0.5,5.0,1.0},{-0.5,5.0,8.0},{-0.5,-5.0,8.0},{0.5,-5.0,1.0}};
    f2.NumVertices = 4;
    f2.barycentre = {0.0,0.0,4.5};

    f3.id = id2;
    f3.vertices = {{3.0, 3.0, 5.0}, {-3.0, 3.0, 4.0}, {-3.0, -3.0, 4.0}, {3.0, -3.0, 5.0}};
    f3.NumVertices = 4;
    f3.barycentre = {0.0,0.0,4.5};

    mesh.NumFractures = 2;
    mesh.FractureId = {0,1};
    mesh.MapFractures.push_back(f1);
    mesh.MapFractures.push_back(f2);

    findIntersections(mesh);

    EXPECT_EQ(mesh.MapTrace.size(),1);
    EXPECT_EQ(mesh.MapFractures[id0].listPas.size(),1);
    EXPECT_EQ(mesh.MapFractures[id1].listNonpas.size(),1);


}

TEST(FRACTURETEST, Testfindintersections2){
    FractureMesh mesh;
    Fracture f1;
    Fracture f2;
    //Fracture f3;
    unsigned int id0 = 0;
    unsigned int id1 = 1;
    //unsigned int id2 = 2;

    f1.id = id0;
    f1.vertices = {{2.0, 1.0, 0.0}, {-2.0, 1.0, 0.0}, {-2.0, -2.0, 0.0}, {2.0, -2.0, 0.0}};
    f1.NumVertices = 4;
    f1.barycentre = {0.0,-0.5,0.0};

    f2.id = id1;
    f2.vertices = {{1.0,-3.0,0.0},{-1.0, -3.0, 4.0},{-1.0, 2.0, 4.0},{1.0, 2.0, 0.0}};
    f2.NumVertices = 4;
    f2.barycentre = {0.0,-0.5,2.0};

    mesh.NumFractures = 2;
    mesh.FractureId = {0,1};
    mesh.MapFractures.push_back(f1);
    mesh.MapFractures.push_back(f2);

    findIntersections(mesh);

    EXPECT_TRUE(mesh.MapFractures[1].onEdge[1]);

}


//TEST funzione void defNewTrace(Trace& t, const double& d1, const double& d2, Fracture& f1, Fracture& f2, FractureMesh& fm)
//CASO TRACCIA NON PASSANTE
TEST(FRACTURETEST, defNewTrace1){
    //uso piani di coordinate f1 = {{1.0, 2.0, 3.0}, {-1.0, 0.0, 2.0}, {-1.0, 0.0, 7.0}, {1.0, 2.0, 7.0}};
    // e f2 = {{0,2,1},{0,-2,1},{0,-2,6},{0,2,6}}
    // i punti di intersezione tra piano e retta sono intersectionf1 = {(0,1,2.5),(0,1,7)} e intersectionf2 = {(0,1,6),(0,1,1)}
    Trace t;
    Fracture f1;
    Fracture f2;
    FractureMesh fm;
    double d1 = 4.5; // distanza punti di intersezione f1
    double d2 = 5.0; // distanza punti di intersezione f2

    t.coordTrace = {{0.0,1.0,2.5},{0.0,1.0,6.0}}; // coordinate traccia
    t.len = dist(t.coordTrace[0],t.coordTrace[1]);; // 3.5

    defNewTrace(t,d1,d2,f1,f2,fm);

    EXPECT_EQ(f1.listPas.size(),0);
    EXPECT_EQ(f2.listPas.size(),0);
    EXPECT_EQ(f1.listNonpas.size(),1);
    EXPECT_EQ(f2.listNonpas.size(),1);
}

//CASO TRACCIA PASSANTE (per f1)
TEST(FRACTURETEST, defNewTrace2){
    //uso piani di coordinate f1 = {{1.0, 2.0, 3.0}, {-1.0, 0.0, 2.0}, {-1.0, 0.0, 7.0}, {1.0, 2.0, 7.0}};
    // e f2 = {{0,2,1},{0,-2,1},{0,-2,8},{0,2,8}}
    // i punti di intersezione tra piano e retta sono intersectionf1 = {(0,1,2.5),(0,1,7)} e intersectionf2 = {(0,1,8),(0,1,1)}
    Trace t;
    Fracture f1;
    Fracture f2;
    FractureMesh fm;
    double d1 = 4.5; // distanza punti di intersezione f1
    double d2 = 7.0; // distanza punti di intersezione f2

    t.coordTrace = {{0.0,1.0,2.5},{0.0,1.0,7.0}}; // coordinate traccia
    t.len = dist(t.coordTrace[0],t.coordTrace[1]);; // 4.5

    defNewTrace(t,d1,d2,f1,f2,fm);

    EXPECT_EQ(f2.listNonpas.size(),1);
    EXPECT_EQ(f1.listPas.size(),1);
    EXPECT_EQ(f1.listNonpas.size(),0);
    EXPECT_EQ(f2.listPas.size(),0);

}





/*
TEST(FRACTURETEST, TestPlotParaviewTraces){

    MatrixXd points = MatrixXd::Zero(3, 4);

    points << 0.0, 1.0, 1.0, 0.0,
        0.0, 0.0, 1.0, 1.0,
        0.0, 0.0, 0.0, 0.0;

    const std::vector<std::vector<unsigned int>> triangles =
        {
            { 0, 1, 3 },
            { 1, 2, 3 }
        };

    Eigen::VectorXi materials(2);
    materials << 0, 1;

    Gedim::UCDUtilities exporter;

    const std::vector<Gedim::UCDProperty<double> > points_properties = {};
    const std::vector<Gedim::UCDProperty<double> > polygons_properties = {};

    exporter.ExportPolygons("./Geometry2Ds.inp",
                            points,
                            triangles,
                            points_properties,
                            polygons_properties,
                            materials);


}
//********************************
TEST(FRACTURETEST, TestPlotParaviewfractures){

    MatrixXd points = MatrixXd::Zero(3, 4);

    points << 0.0, 1.0, 1.0, 0.0,
        0.0, 0.0, 1.0, 1.0,
        0.0, 0.0, 0.0, 0.0;

    vector<vector<unsigned int>> pol_vertices = { {0, 1, 2, 3}};

    Polygons polygons(points,
                      pol_vertices);
    std::vector<std::vector<unsigned int>> triangles;
    Eigen::VectorXi materials;

    polygons.GedimInterface(triangles, materials);

    Gedim::UCDUtilities exporter;

    exporter.ExportPolygons("./Geometry2Ds.inp",
                            points,
                            triangles,
                            {},
                            {},
                            materials);
}
//********************************


*/
}

#endif // TEST_HPP
