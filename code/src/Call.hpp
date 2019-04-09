#ifndef CALL_H
#define CALL_H

#include "option.hpp"

class Call: public Option{
    
    public:
        Call(double _strike, double _horizon): Option(_strike, _horizon) {};
        double boundary_left(double t, double x, double rate);
        double boundary_right(double t, double x, double rate);
        double payoff(double x);
};

#endif