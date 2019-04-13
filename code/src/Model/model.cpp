#include "model.hpp"

Model::Model(std::function<double (double,double)> _rate,
             std::function<double (double,double)> _sigma,
             double _S0)
    : rate(_rate), sigma(_sigma), S0(_S0){};

double Model::S(){ return S0;}
