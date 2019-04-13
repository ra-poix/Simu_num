#ifndef BINARY_H
#define BINARY_H

#include "Option.hpp"

class Binary: public Option{

    private: 
        double strike;
        double horizon;
        
    public:
        Binary(double _strike, double _horizon): strike(_strike), horizon(_horizon){};
        ~Binary();
        double boundary_left(double t, double x, double rate);
        double boundary_right(double t, double x, double rate);
        double payoff(double x);

        double Strike(){ return strike; };
        double Horizon(){ return horizon; };
};

#endif