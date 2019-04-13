//
//  euler.cpp
//  
//
//  Created by Côme RIMBERT on 13/04/2019.
//
/*
#include "euler.hpp"
#include <cmath>

double euler::solution(int N){
    double Y=x.S();
    double T=z.Horizon();
    double pas=T/N;
    double time=0;
    while(time<T){
        Y+= Y + x.rate(time,Y)*pas + x.sigma(time,Y)*gen.normale()*sqrt(pas);
        time+=pas;
    }
    return Y;
}
*/