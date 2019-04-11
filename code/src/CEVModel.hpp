#ifndef CEVMODEL_H
#define CEVMODEL_H
#include "model.hpp"
#include <algorithm>
#include <math.h>

class CEVModel : public Model{

    public:
        CEVModel(double sigma, double r, double s, double NU): sigma(sigma), r(r), S0(s), nu(NU){}
        ~CEVModel(){};

        double S(){ return S0; };
        double Sigma(double t, double x);
        double Rate(double t, double x);

    private:
        double sigma;
        double r;
        double S0;
        double nu;

};
#endif