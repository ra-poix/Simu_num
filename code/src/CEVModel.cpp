#include "CEVModel.hpp"

double CEVModel::Sigma(double t, double x){ 
    return pow(x,nu) / sqrt(1+x*x);}
double CEVModel::Rate(double t, double x){ 
    return r;
}