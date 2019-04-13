#ifndef OPTION_H
#define OPTION_H

class Option{

    public:
        virtual double boundary_left(double t, double x, double rate){ return 0;};
        virtual double boundary_right(double t, double x, double rate){ return 0;};
        virtual double payoff(double x){ return x; };

        virtual double Horizon () const { return horizon;};
        
    private:
        double horizon;
};

#endif