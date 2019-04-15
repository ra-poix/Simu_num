#include "Model/BSmodel.hpp"
#include "Model/euler.hpp"
#include "Model/CEVModel.hpp"

#include "Option/Call.hpp"
#include "Option/Binary.hpp"

#include "MeanVar.hpp"
#include "utils.hpp"

#include "Generateur/Generateur.hpp"
#include "Generateur/RandomGenerator.hpp"
#include "Generateur/LowDiscrepancyGenerator.hpp"
#include "Generateur/LowDiscrepancySequence.hpp"

#include <cstdio>
#include <ctime>
#include <iostream>
#include <fstream>
#include <new>
#include <string.h>

//pointeur sur fonction
//Utilisation : function mon_pointeur[] = {maFonction1, maFonction 2, maFonction3, ...}
//type de retour : double
//params d'entrée : Option, Model, double, double
typedef double (*functions) (Option, Model, double*, double);

int MAX_SIZE_MC = 10000000;
int MAX_SIZE_MC_EULER = 100000;
int EULER_PRECISION = 1000;

bool enough_precise(MeanVar *mv, int size_mv) {
    for(int i = 0 ; i < size_mv; i++){
        if(!mv[i].is_enough_precise())
            return false;
    }
    return true;
}
/*
int MonteCarlo(Option &o, Model &m, 
                Generateur &g, 
                functions func, 
                std::function<double (int)> pas,
                MeanVar *mv, bool is_euler){
    int compteur = 0;
    double *X;
    while(compteur < MAX_SIZE_MC){
        compteur ++;
        if(compteur % 100000 == 0 && mv -> is_enough_precise())
            break;
        else if (compteur%1000000 ==0)
            std::cout << " MC running : " << compteur << " simulations" << std::endl;
        if(is_euler)
            *X = g.normale();
        else{
            euler scheme(o,m,g);
            X = scheme.solution(EULER_PRECISION,true, 0.01);
        }
      //  std::cout << pas(compteur) << " " << std::endl;fflush(stdout);
        mv -> maj( func(o,m,X,pas(compteur)) );
  //      mv -> maj( func(o,m,-X,pas(compteur)) );

    }
    return compteur;
};*/

//MC direct
int MonteCarlo(Option &o, BSmodel &m, 
                Generateur &g, 
                functions func, 
                std::function<double (int)> pas,
                MeanVar *mv){
    int compteur = 0;
    while(compteur < MAX_SIZE_MC){
        compteur ++;
        if(compteur % 100000 == 0 && mv -> is_enough_precise())
            break;
        else if (compteur%1000000 ==0)
            std::cout << " MC running : " << compteur << " simulations" << std::endl;
        double X = g.normale();
      //  std::cout << pas(compteur) << " " << std::endl;fflush(stdout);
        mv -> maj( func(o,m,&X,pas(compteur)) );
  //      mv -> maj( func(o,m,-X,pas(compteur)) );

    }
    return compteur;
};


//MonteCarlo pour schema d'euler
int MonteCarlo(Option &o, CEVModel &m, 
                Generateur &g, 
                functions func,
                std::function<double (int)> pas,
                MeanVar *mv){
    int compteur = 0;
    while(compteur < MAX_SIZE_MC_EULER){
        compteur ++;
        if(compteur % 1000== 0 && mv -> is_enough_precise())
            break;
        else if (compteur%100000 ==0)
            std::cout << " MC running : " << compteur << " simulations " << std::endl;
        euler scheme(o,m,g);
        double *X = scheme.solution(EULER_PRECISION,true, 0.01);
        // printf("MonteCarlo g = %f \n",g);
        // std::cout << *X << std::endl;fflush(stdout);
         //std::cout << " e " << func(o,m,X,pas(compteur)) ;
        mv -> maj( func(o,m,X,pas(compteur)) );        
    }
    return compteur;
};

/*void run(std::ofstream os, 
        CEVModel b, Option c, functions f,Generateur g,MeanVar *mv,
        double precision,double S,std::function<double (int)> pas,){
    os << "S" << ";" << "mean" << ";" << "var"<< ";" << "nb_sim" << ";" << "time" << std::endl;
    std::cout << "\n START \n";
    for(S = 20 ; S < 180 ; S++){
        std::clock_t start;
        double duration;
        start = std::clock();
        int compteur = MonteCarlo(c, , g, f, pas, mv);
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        std::cout << "rate: " << b.rate(0.0,0.0) << ", vol: " << b.sigma(0.0,0.0) << ", strike: "<< c.Strike() << ", horizon: " << c.Horizon() << std::endl;
        std::cout<< "S = " << S << ", mean = "<< mv->mean() << ", var = " << mv -> var() <<", compteur = " << compteur << ", time = " << duration <<  std::endl;
        os << S << ";" << mv->mean() << ";" << mv->var() << ";" << compteur << ";" << duration << std::endl;
        *mv = MeanVar(precision);
    }
    delete mv;
    os.close();
}*/

int main(){

    std::ofstream aux;
    double S;
    double r = 0.02;
    double vol = 0.2;
    double nu = 0.2;
   // BSmodel b(r, vol, S); 
    CEVModel b(r,S,nu);

    bool is_euler = true;

    double K = 100;
    double T = 5;
    Call c(K,T);

    /*Halton h(2); 
    Generateur g = LowDiscrepancyGenerator(h);*/

    Generateur g = RandomGenerator();

   std::function<double (int)> pas_const =  [&] (int i){ return 0.001;};
   std::function<double (int)> pas_dcr = [&] (int i){ 
       return 1 / pow(i,0.25);
    };

    functions f  = FD_gamma_cev;

    double precision = 0.001;
    MeanVar *mv = new MeanVar;
    *mv = MeanVar(precision);

    char *output_name = "FD_gamma_cev_mat5.csv";
    std::ofstream fdm_out = std::ofstream(output_name);
    fdm_out << "S" << ";" << "mean" << ";" << "var"<< ";" << "nb_sim" << ";" << "time" << std::endl;

    std::cout << "\n START \n";

    for(S = 20 ; S < 180 ; S++){

        std::clock_t start;
        double duration;
        start = std::clock();

        //BSmodel b(r,vol,S);
        CEVModel b(r,S,nu);

        int compteur = MonteCarlo(c, b, g, f, pas_dcr, mv);
        //        int compteur = MonteCarlo(c, b, g, f, pas_dcr, mv, is_euler);
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

        std::cout << "rate: " << b.rate(0.0,0.0) << ", vol: " << b.sigma(0.0,0.0) << ", strike: "<< c.Strike() << ", horizon: " << c.Horizon() << std::endl;
        std::cout<< "S = " << S << ", mean = "<< mv->mean() << ", var = " << mv -> var() <<", compteur = " << compteur << ", time = " << duration <<  std::endl;
        fdm_out << S << ";" << mv->mean() << ";" << mv->var() << ";" << compteur << ";" << duration << std::endl;
        *mv = MeanVar(precision);
    }

    delete mv;

    //fdm_out.close();

    std::cout << "\nFIN MC\n";
    return 0;
}