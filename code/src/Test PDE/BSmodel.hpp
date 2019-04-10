#ifndef BSMODEL_H
#define BSMODEL_H
#include "model.hpp"

class BSmodel : public Model{

    public:
        BSmodel(double sigma, double r, double s): sigma(sigma), r(r), S0(s){}
        ~BSmodel(){};

        double S(){ return S0; };
        double Sigma(double t, double x);
        double Rate(double t, double x);

    private:
        double sigma;
        double r;
        double S0;

};
#endif