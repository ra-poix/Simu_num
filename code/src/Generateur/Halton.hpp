#ifndef HALTON_H
#define HALTON_H

#include "LowDiscrepancySequence.hpp"
#include <iostream>

class Halton: public LowDiscrepancySequence{
    public:
        Halton(long _dim, long index = 1);
        ~Halton(){};

        std::vector<double> next();

    private:
        std::vector<long> p;
        long index;

        void set_params();
        void find_primes(long _dim);
        std::vector<long> decompose(long n, long p);
        bool is_prime(long p);
        double phi(long p);

};



#endif