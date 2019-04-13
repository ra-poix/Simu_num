#ifndef CEVMODEL_H
#define CEVMODEL_H
#include "model.hpp"
#include <algorithm>
#include <math.h>

class CEVModel : public Model{
    
    typedef double (*function) (double,double);


    public:
        CEVModel(double r, double s, double NU);
        ~CEVModel(){};
 
    private:
        double S0;
        double nu;

};
#endif