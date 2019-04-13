//
//  euler.hpp
//  
//
//  Created by Côme RIMBERT on 13/04/2019.
//

#ifndef euler_hpp
#define euler_hpp

#include "model.hpp"
#include "Generateur.hpp"
#include "option.hpp"

class euler{
public:
    //euler(Option z,Model x,Generateur gen): z(z) x(x), gen(gen) {};
    ~euler(){};
    double solution(int N);
    
private:
    Model x;
    Generateur gen;
    Option z;
};
#endif /* euler_hpp */
