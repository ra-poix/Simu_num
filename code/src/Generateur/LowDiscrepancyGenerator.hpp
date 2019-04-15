#ifndef HALTONGENERATOR_H
#define HALTONGENERATOR_H

#include "Generateur.hpp"
#include "Halton.hpp"
#include "LowDiscrepancySequence.hpp"
#include <cmath>

class LowDiscrepancyGenerator: public Generateur{

    public:
        LowDiscrepancyGenerator(LowDiscrepancySequence _seq):
            Generateur( [this] () { 
                std::vector<double> a = seq.next();
                return sqrt(-2 * log (a[0]) ) * cos(2*std::atan(1)*4*a[1]);
            }),
            seq(_seq) {};

    private:
        LowDiscrepancySequence seq;

};
#endif