#ifndef UTILS_HPP
#define UTILS_HPP

#include "dfn.hpp"

using namespace std;
using namespace Eigen;

namespace FractureLibrary
{

bool ImportFR_data(const string &filename,
                   FractureMesh& mesh);

}
#endif // UTILS_HPP
