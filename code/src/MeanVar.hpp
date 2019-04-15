#ifndef MEANVAR_H
#define MEANVAR_H

#include <math.h>

class MeanVar{

    private:
        double precision;
        int n;
        double x;
        double xx;
    
    public:
        MeanVar(){};
        MeanVar(double _precision):precision(_precision),n(0), x(0), xx(0) {};
 
        void maj(double y){
            n++;
            x+=y;
            xx+=pow(y,2);
        };

        double mean() const {
            return x/n;
        };

        double var() const{
            return xx/n - pow(mean(),2);
        }

        bool is_enough_precise() const { 
            double aux = 1.96 * sqrt( var() /n ) ;
            return ( 1.96 * sqrt(var() /n) < precision); 
        }

        int size() const { return n;};
        int sum_x() const { return x;};
        int sum_xx() const { return xx;};
};

#endif