#ifndef UCD_test_HPP
#define UCD_test_HPP

#include "UCDUtilities.hpp"
#include <gtest/gtest.h>

// *************************
// ESPORTO SU PARAVIEW FRATTURA 1 DI FR3
TEST(TestParaviewfr_3, UCDUtilities_Test0Ds_FR3_1)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 8)<< 0.0, 1.0, 1.0, 0.0, 0.8,0.8,0.8,0.0,
                                    0.0, 0.0, 1.0, 1.0,0.0,1.0,0.5,0.5,
                                    0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0).finished();

    exporter.ExportPoints(exportFolder + "/Geometry0Ds_fr3_1.inp",
                          points);

}
// *************************
TEST(TestParaviewfr_3, UCDUtilities_Test1Ds_FR3_1)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 8)<< 0.0, 1.0, 1.0, 0.0, 0.8,0.8,0.8,0.0,
                                    0.0, 0.0, 1.0, 1.0,0.0,1.0,0.5,0.5,
                                    0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0).finished();
    const Eigen::MatrixXi edges = (Eigen::MatrixXi(2, 10)<< 1, 0, 4, 2, 5,4,6,3,7,7,
                                   2,4,1,5,3,6,5,7,0,6).finished();

    exporter.ExportSegments(exportFolder + "/Geometry1Ds_fr3_1.inp",
                            points,
                            edges);
}
// *************************
TEST(TestParaviewfr_3, UCDUtilities_Test2Ds_FR3_1)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 8)<< 0.0, 1.0, 1.0, 0.0, 0.8,0.8,0.8,0.0,
                                    0.0, 0.0, 1.0, 1.0,0.0,1.0,0.5,0.5,
                                    0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0).finished();
    const std::vector<std::vector<unsigned int>> polygons =
        {
            { 0, 4, 6,7 },
            { 6, 5, 3, 7 },
            { 4, 1, 2 },
            { 4, 2, 5 },
            { 4, 5, 6 }
        };

    exporter.ExportPolygons(exportFolder + "/Geometry2Ds_fr3_1.inp",
                            points,
                            polygons);
}
// *************************



// ESPORTO SU PARAVIEW FRATTURA 2 DI FR3
TEST(TestParaviewfr_3, UCDUtilities_Test0Ds_FR3_2)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 6)<< 0.8, 0.8, 0.8, 0.8, 0.8,0.8,
                                    0.0, 0.0, 1.0, 1.0,0.0, 1.0,
                                    -0.1,0.29999,0.29999,-0.1,0.0,0.0).finished();

    exporter.ExportPoints(exportFolder + "/Geometry0Ds_fr3_2.inp",
                          points);

}
// *************************
TEST(TestParaviewfr_3, UCDUtilities_Test1Ds_FR3_2)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 6)<< 0.8, 0.8, 0.8, 0.8, 0.8,0.8,
                                    0.0, 0.0, 1.0, 1.0,0.0, 1.0,
                                    -0.1,0.29999,0.29999,-0.1,0.0,0.0).finished();

    const Eigen::MatrixXi edges = (Eigen::MatrixXi(2, 7)<< 1,3,0,4,2,5,5,
                                   2,0,4,1,5,3,4).finished();

    exporter.ExportSegments(exportFolder + "/Geometry1Ds_fr3_2.inp",
                            points,
                            edges);
}
// *************************
TEST(TestParaviewfr_3, UCDUtilities_Test2Ds_FR3_2)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 6)<< 0.8, 0.8, 0.8, 0.8, 0.8,0.8,
                                    0.0, 0.0, 1.0, 1.0,0.0, 1.0,
                                    -0.1,0.29999,0.29999,-0.1,0.0,0.0).finished();

    const std::vector<std::vector<unsigned int>> polygons =
        {
            {0,4,5,3},
            {4,1,2,5}
        };

    exporter.ExportPolygons(exportFolder + "/Geometry2Ds_fr3_2.inp",
                            points,
                            polygons);
}
// *************************



// ESPORTO SU PARAVIEW FRATTURA 3 DI FR3
TEST(TestParaviewfr_3, UCDUtilities_Test0Ds_FR3_3)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 6)<< -0.237777,0.316183700,0.316183700,-0.23777,0.316183700,-0.23777,
                                    0.5,0.5,0.5,0.5,0.5,0.5,
                                    -0.344440, -0.344440,0.452838899,0.452838899,0.0,0.0).finished();

    exporter.ExportPoints(exportFolder + "/Geometry0Ds_fr3_3.inp",
                          points);

}
// *************************
TEST(TestParaviewfr_3, UCDUtilities_Test1Ds_FR3_3)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 6)<< -0.237777,0.316183700,0.316183700,-0.23777,0.316183700,-0.23777,
                                    0.5,0.5,0.5,0.5,0.5,0.5,
                                    -0.344440, -0.344440,0.452838899,0.452838899,0.0,0.0).finished();

    const Eigen::MatrixXi edges = (Eigen::MatrixXi(2,7)<< 0,2,1,4,3,5,5,
                                   1,3,4,2,5,0,4).finished();

    exporter.ExportSegments(exportFolder + "/Geometry1Ds_fr3_3.inp",
                            points,
                            edges);
}
// *************************
TEST(TestParaviewfr_3, UCDUtilities_Test2Ds_FR3_3)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 6)<< -0.237777,0.316183700,0.316183700,-0.23777,0.316183700,-0.23777,
                                    0.5,0.5,0.5,0.5,0.5,0.5,
                                    -0.344440, -0.344440,0.452838899,0.452838899,0.0,0.0).finished();

    const std::vector<std::vector<unsigned int>> polygons =
        {
            {0,1,4,5},
            {4,2,3,5}
        };

    exporter.ExportPolygons(exportFolder + "/Geometry2Ds_fr3_3.inp",
                            points,
                            polygons);
}

// *************************

TEST(TestParaviewfr_10, UCDUtilities_Test0Ds_FR10)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 24)<< 0.679497,0.209594,0.0770272,0.54693,0.422841,0.554138,0.160363,0.492678,0.436999,0.594717,0.407083,0.329511, 0.364069,0.446922,0.450525,0.247976,0.241961,0.235293,0.431609,0.433893,0.240595,0.479352,0.4272,0.352669,
                                    0.515669,0.993894,0.863634,0.385409,0.77687,0.392492,0.945519,0.440621,0.735424,0.432365,0.792907,0.606678,0.636022,0.706373,0.709432,0.689657,0.821544,0.967739,0.74578,0.744516,0.851491,0.719358,0.507259,0.626342,
                                    0.190545,0.190545,0.597042,0.597042,0.190545,0.574938,0.341506,0.597042,0.231993,0.450509,0.190545,0.597042,0.498147,0.261045,0.250734,0.597042,0.404251,0.190545,0.224526,0.2229,0.360475,0.190545,0.597042,0.530771).finished();

    exporter.ExportPoints(exportFolder + "/Geometry0Ds_fr10.inp",
                          points);

}
// *************************
TEST(TestParaviewfr_10, UCDUtilities_Test1Ds_FR10)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;

    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 24)<< 0.679497,0.209594,0.0770272,0.54693,0.422841,0.554138,0.160363,0.492678,0.436999,0.594717,0.407083,0.329511, 0.364069,0.446922,0.450525,0.247976,0.241961,0.235293,0.431609,0.433893,0.240595,0.479352,0.4272,0.352669,
                                    0.515669,0.993894,0.863634,0.385409,0.77687,0.392492,0.945519,0.440621,0.735424,0.432365,0.792907,0.606678,0.636022,0.706373,0.709432,0.689657,0.821544,0.967739,0.74578,0.744516,0.851491,0.719358,0.507259,0.626342,
                                    0.190545,0.190545,0.597042,0.597042,0.190545,0.574938,0.341506,0.597042,0.231993,0.450509,0.190545,0.597042,0.498147,0.261045,0.250734,0.597042,0.404251,0.190545,0.224526,0.2229,0.360475,0.190545,0.597042,0.530771).finished();


    const Eigen::MatrixXi edges = (Eigen::MatrixXi(2, 37)<< 3,1,6,7,5,9,4,7,5,13,12,9,14,13,2,15,12,16,16,10,17,10,18,8,19,19,17,20,18,0,21,19,11,22,12,23,23
                                   ,5,6,2,3,9,0,10,12,13,8,13,14,8,14,15,11,16,6,15,17,1,18,8,19,4,18,20,16,20,21,4,21,22,7,23,11,22).finished();

    exporter.ExportSegments(exportFolder + "/Geometry1Ds_fr10.inp",
                            points,
                            edges);
}
// *************************
TEST(TestParaviewfr_10, UCDUtilities_Test2Ds_FR10)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 24)<< 0.679497,0.209594,0.0770272,0.54693,0.422841,0.554138,0.160363,0.492678,0.436999,0.594717,0.407083,0.329511, 0.364069,0.446922,0.450525,0.247976,0.241961,0.235293,0.431609,0.433893,0.240595,0.479352,0.4272,0.352669,
                                    0.515669,0.993894,0.863634,0.385409,0.77687,0.392492,0.945519,0.440621,0.735424,0.432365,0.792907,0.606678,0.636022,0.706373,0.709432,0.689657,0.821544,0.967739,0.74578,0.744516,0.851491,0.719358,0.507259,0.626342,
                                    0.190545,0.190545,0.597042,0.597042,0.190545,0.574938,0.341506,0.597042,0.231993,0.450509,0.190545,0.597042,0.498147,0.261045,0.250734,0.597042,0.404251,0.190545,0.224526,0.2229,0.360475,0.190545,0.597042,0.530771).finished();

    const std::vector<std::vector<unsigned int>> polygons =
        {
         /*{12,7,3,5,13},
            {8,13,14 },
            {13,5,9,14},
            {6,2,15,16},
            {4,10,18,19},
        {18,8,19},
        {10,17,20,18},
        {20,16,12,13,8,18},
        {0,21,19,8,14,9},
        {21,4,19},
        {17,1,6,16,20},
        {11,22,23},
        {22,7,12,23},
        {15,11,23,12,16}
        };*/
         {12, 7, 3}, {12, 3, 5}, {12, 5, 13}, // {12, 7, 3, 5, 13}
         {8, 13, 14}, // {8, 13, 14}
         {13, 5, 9}, {13, 9, 14}, // {13, 5, 9, 14}
         {6, 2, 15}, {6, 15, 16}, // {6, 2, 15, 16}
         {4, 10, 18}, {4, 18, 19}, // {4, 10, 18, 19}
         {18, 8, 19}, // {18, 8, 19}
         {10, 17, 20}, {10, 20, 18}, // {10, 17, 20, 18}
         {20, 16, 12}, {20, 12, 13}, {20, 13, 8}, {20, 8, 18}, // {20, 16, 12, 13, 8, 18}
         {0, 21, 19}, {0, 19, 8}, {0, 8, 14}, {0, 14, 9}, // {0, 21, 19, 8, 14, 9}
         {21, 4, 19}, // {21, 4, 19}
         {17, 1, 6}, {17, 6, 16}, {17, 16, 20}, // {17, 1, 6, 16, 20}
         {11, 22, 23}, // {11, 22, 23}
         {22, 7, 12}, {22, 12, 23}, // {22, 7, 12, 23}
         {15, 11, 23}, {15, 23, 12}, {15, 12, 16} };

    exporter.ExportPolygons(exportFolder + "/Geometry2Ds_fr10.inp",
                            points,
                            polygons);
}
#endif // UCD_test_HPP
