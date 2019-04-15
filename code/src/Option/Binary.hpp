#ifndef BINARY_H
#define BINARY_H

#include "Option.hpp"

class Binary: public Option{
        
    public:
        Binary(double _strike, double _horizon);
        ~Binary(){};
        double boundary_left(double t, double x, double rate);
        double boundary_right(double t, double x, double rate);
        double Horizon() const { return horizon;};
        double Strike() const { return strike;}
    private:
        double horizon;
        double strike;

};

#endif