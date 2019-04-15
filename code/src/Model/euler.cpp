//
//  euler.cpp
//  
//
//  Created by Côme RIMBERT on 13/04/2019.
//

#include "euler.hpp"

#include <fstream>


double* euler::solution(int N, bool is_fd, double delta){
    double *result;
    if(is_fd)
        result = new double[1];
    else
        result = new double[3];
    
    double Y1= x.S();
    double Y2= Y1 + delta;
    double Y3= Y1 - delta;

    double T=z.Horizon();
    double pas=T/N;
    double time=0;

    while(time<T){

        double norm = gen.normale();
        Y1 += Y1 * ( x.rate(time,Y1)*pas + x.sigma(time,Y1)*norm*sqrt(pas) );

        if(is_fd){
            Y2 += Y2 * ( x.rate(time,Y2)*pas + x.sigma(time,Y2)*norm*sqrt(pas) );
            Y3 += Y3 * ( x.rate(time,Y3)*pas + x.sigma(time,Y3)*norm*sqrt(pas) );
        }
        time+=pas;   
    }
    result[0] = Y1;
    if(is_fd){
        result[1] = Y2;
        result[2] = Y3;
    }
    return result;
}

