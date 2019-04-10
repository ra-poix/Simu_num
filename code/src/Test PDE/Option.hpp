#ifndef OPTION_H
#define OPTION_H

class Option{

    public:
        virtual double boundary_left(double t, double x, double rate) = 0;
        virtual double boundary_right(double t, double x, double rate) = 0;
        virtual double payoff(double x) = 0;
        
};

#endif