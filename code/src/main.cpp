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

#include <cmath>
#include <iostream>
#include <fstream>
#include <new>

//pointeur sur fonction
//Utilisation : function mon_pointeur[] = {maFonction1, maFonction 2, maFonction3, ...}
//type de retour : double
//params d'entrée : Option, Model, double, double
typedef double (*functions) (Option, Model, double, double);

bool enough_precise(MeanVar *mv, int size_mv) {
    for(int i = 0 ; i < size_mv; i++){
        if(!mv[i].is_enough_precise())
            return false;
    }
    return true;
}

//MC direct
void MonteCarlo(Option &o, BSmodel &m, Generateur &g, functions func, MeanVar *mv){
    int compteur = 0;
    while(compteur < 1000000){
        compteur ++;
        if(compteur % 100000 == 0 && mv -> is_enough_precise())
            break;
        //std::cout << "oijojioj" ;
        double X = g.normale();
       // std::cout << X << " " << std::endl;fflush(stdout);
        mv -> maj( func(o,m,X,1) );
        mv -> maj( func(o,m,-X,1) );

    }

    std::cout << "compteur = " << compteur << std::endl;
};

//MonteCarlo pour schema d'euler
void MonteCarlo(Option &o, CEVModel &m, RandomGenerator &g, functions func, MeanVar *mv){
    int compteur = 0;
    while(compteur < 1){
        compteur ++;
        if(compteur % 100000 == 0 && mv -> is_enough_precise())
            break;
        euler scheme(o,m,g);
        double X=scheme.solution(1000);
        // printf("MonteCarlo g = %f \n",g);
        // std::cout << X << std::endl;fflush(stdout);
        mv -> maj( func(o,m,X,1) );        
    }
    
    std::cout << "compteur = " << compteur << std::endl;
};

int main(){

    functions f  = Tangent_delta_bs;
    int size_f = 1;

    double S;
    double r = 0.02;
    double vol = 0.4;
    BSmodel b(r, vol, S);

    double nu = 0.5;
    //CEVModel b(r,S,nu);

    double K = 100;
    double T = 1;
    Call c(K,T);

    double precision = 0.001;
    MeanVar *mv = new MeanVar;
    *mv = MeanVar(precision);
    std::cout << "oazidja" << std::endl;

    Halton h(2);
    Generateur g = LowDiscrepancyGenerator(h);

   // Generateur g = RandomGenerator();

    std::ofstream fdm_out("fdm.csv");
    fdm_out << "S" << ";" << "Y" << std::endl;
    for(S = 100 ; S < 101 ; S++){
        BSmodel b(r,vol,S);
        //CEVModel b(r,S,nu);
        MonteCarlo(c, b, g, f, mv);
        std::cout<< "S = " << S << " mean = "<< mv->mean() << ", var = " << mv -> var() << std::endl;
        fdm_out << S << ";" << mv->mean() << std::endl;
        *mv = MeanVar(precision);
    }

    delete mv;

    std::cout << "\nFIN MC\n";
    std::cout << "mean = "<< mv[0].mean() << ", var = " << mv[0].var() << std::endl;

    return 0;
}