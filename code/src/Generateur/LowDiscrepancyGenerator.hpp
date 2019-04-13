#ifndef HALTONGENERATOR_H
#define HALTONGENERATOR_H

#include "Generateur.hpp"
#include "Halton.hpp"
#include "LowDiscrepancySequence.hpp"
#include <cmath>

class LowDiscrepancyGenerator: public Generateur{

    public:
        LowDiscrepancyGenerator(LowDiscrepancySequence _seq = Halton(2)): Generateur(), seq(_seq) {};
        double normale(){
            std::vector<double> a = seq.next();
            return sqrt(-2 * log (a[0]) ) * cos(2*std::atan(1)*4*a[1]);
        };

    private:
        LowDiscrepancySequence seq;

};
#endif