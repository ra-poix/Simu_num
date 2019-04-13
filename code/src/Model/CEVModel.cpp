#include "CEVModel.hpp"

CEVModel::CEVModel(double r, double s, double NU) :
    Model(  [r] (double t, double x) { return r;}, 
            [NU] (double t, double x) { return pow(x,NU) / sqrt(1 + pow(x,2)) ;}, 
            s) {}
