#ifndef CALL_H
#define CALL_H

#include "option.hpp"
#include <math.h>

class Call: public Option{

    private: 
        double strike;
        double horizon;
        
    public:
        Call(double _strike, double _horizon);

        ~Call(){};
        double boundary_left(double t, double x, double rate);
        double boundary_right(double t, double x, double rate);
      //  double payoff(double x);

        double Strike(){ return strike; };
        double Horizon(){ return horizon; };
};

#endif