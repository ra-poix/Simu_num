//
//  euler.cpp
//  
//
//  Created by Côme RIMBERT on 13/04/2019.
//

#include "euler.hpp"

#include <fstream>


double euler::solution(int N){
    double Y=x.S();
    double T=z.Horizon();
    double pas=T/N;
    double time=0;

    std::ofstream e("euler.csv");
    e << "T" << ";" << "Y" << std::endl;
    e << time << ";" << Y << std::endl;
    while(time<T){
        Y+= Y* ( x.rate(time,Y)*pas + x.sigma(time,Y)*gen.normale()*sqrt(pas) );
        e << time << ";" << Y << std::endl;
        time+=pas;   
    }
    return Y;
    e.close();
}

