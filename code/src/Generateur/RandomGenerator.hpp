#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include "Generateur.hpp"
#include <random>
#include <ctime>

class RandomGenerator : public Generateur{

    public:
        RandomGenerator(double seed = std::time(nullptr)):
            Generateur(  [this] () { return distribution(gen);} ),
            gen(seed),
            distribution(){} ;

        
    private:
        std::mt19937_64 gen;
        std::normal_distribution<double> distribution;

};

#endif