#ifndef UCD_test_HPP
#define UCD_test_HPP

#include "UCDUtilities.hpp"
#include <gtest/gtest.h>

// *************************
TEST(TestUCDUtilities, UCDUtilities_Test0Ds)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 8)<< 0.0, 1.0, 1.0, 0.0, 0.8,0.8,0.8,0.0,
                                    0.0, 0.0, 1.0, 1.0,0.0,1.0,0.5,0.5,
                                    0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0).finished();

    exporter.ExportPoints(exportFolder + "/Geometry0Ds.inp",
                          points);

}
// *************************
TEST(TestUCDUtilities, UCDUtilities_Test1Ds)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 8)<< 0.0, 1.0, 1.0, 0.0, 0.8,0.8,0.8,0.0,
                                    0.0, 0.0, 1.0, 1.0,0.0,1.0,0.5,0.5,
                                    0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0,0.0).finished();
    const Eigen::MatrixXi edges = (Eigen::MatrixXi(2, 10)<< 1, 0, 4, 2, 5,4,6,3,7,7,
                                   2,4,1,5,3,6,5,7,0,6).finished();

    exporter.ExportSegments(exportFolder + "/Geometry1Ds.inp",
                            points,
                            edges);
}
// *************************
TEST(TestUCDUtilities, UCDUtilities_Test2Ds)
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

    exporter.ExportPolygons(exportFolder + "/Geometry2Ds.inp",
                            points,
                            polygons);
}
// *************************

#endif // UCD_test_HPP
