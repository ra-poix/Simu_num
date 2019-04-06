#ifndef BSMODEL_H
#define BSMODEL_H
#include "model.hpp"

class BSmodel : public Model{

    public:
        BSmodel();
        BSmodel(double sigma, double r): sigma(sigma), r(r){}

        ~BSmodel();

        double Sigma();
        double Rate();

    private:
        double sigma;
        double r;

};
#endif