//
//  euler.hpp
//  
//
//  Created by Côme RIMBERT on 13/04/2019.
//

#ifndef euler_hpp
#define euler_hpp

#include "model.hpp"
#include "../Generateur/Generateur.hpp"
#include "../Generateur/RandomGenerator.hpp"

#include "../Option/option.hpp"

#include <cmath>

class euler{
public:
    euler(Option z,Model x,RandomGenerator &gen): z(z), x(x), gen(gen) {};
    ~euler(){};
    double solution(int N);
    
private:
    Model x;
    RandomGenerator &gen;
    Option z;
};
#endif /* euler_hpp */
