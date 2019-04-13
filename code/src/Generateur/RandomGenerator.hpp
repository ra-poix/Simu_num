#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include "Generateur.hpp"
#include <random>
#include <ctime>

class RandomGenerator : public Generateur{

    public:
        RandomGenerator(double seed = std::time(nullptr)): gen(seed), distribution(){} ;

        double normale(){
            return distribution(gen);
        };

    private:
        std::mt19937_64 gen;
        std::normal_distribution<double> distribution;

};

#endif