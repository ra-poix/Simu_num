#ifndef BSMODEL_H
#define BSMODEL_H
#include "model.hpp"


class BSmodel : public Model{

    public:
        BSmodel(double r,  double vol, double s);
        ~BSmodel(){};
        double S(){ return S0;};

    private:
        double S0;
};

#endif