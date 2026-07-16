#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "random/sfmt.h"
// Third-party RNG (Agner Fog); its type punning trips -Wstrict-aliasing, kept as is
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#include "random/sfmt.cpp"
#pragma GCC diagnostic pop
#include <float.h>
#include <string>
#include "version.h"
using namespace std;
