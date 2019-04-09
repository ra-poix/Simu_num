#ifndef OPTION_H
#define OPTION_H

#include <algorithm>

class Option{

    protected:
        double strike;
        double horizon;

    public:
        Option(double _strike, double _horizon): strike(_strike), horizon(_horizon){};
        virtual double boundary_left(double t, double x, double rate);
        virtual double boundary_right(double t, double x, double rate);
        virtual double payoff(double x);

        double Horizon(){ return strike; };
        double Strike() { return horizon; };
        
};

#endif