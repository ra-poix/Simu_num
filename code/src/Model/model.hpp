#ifndef MODEL_H
#define MODEL_H

#include <functional>
#include <iostream>

class Model{

    public:
    
        std::function<double (double,double)> rate;
        std::function<double (double,double)> sigma;
        Model(){};
        Model(std::function<double (double,double)> _rate, std::function<double (double,double)> _sigma,double _S0);

        virtual double S();

    protected:
        double S0;
};



#endif