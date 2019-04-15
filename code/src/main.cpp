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

#include "PDE/fdm.hpp"

#include <cstdio>
#include <ctime>
#include <iostream>
#include <fstream>
#include <new>
#include <string.h>

typedef double (*functions) (Option, Model, double*, double);

int MAX_SIZE_MC = 10000000; //Pour arrêter le calcul si la variance explose
int MAX_SIZE_MC_EULER = 10000; //Temps de simulation élevé, il faut brider ce param ou la précision du schéma
int EULER_PRECISION = 100; //taille des pas d'un schéma d'euler

int MonteCarlo(Option &o, Model &m, 
                Generateur &g, 
                functions func, 
                std::function<double (int)> pas,
                MeanVar *mv, bool is_euler){
    int compteur = 0;
    double *X;
    int max = (is_euler)?MAX_SIZE_MC_EULER : MAX_SIZE_MC;
    while(compteur < max){
        compteur ++;
        if(compteur % 1000000 == 0 && mv -> is_enough_precise())
            break;
        else if (compteur%1000000 ==0)
            std::cout << " MC running : " << compteur << " simulations" << std::endl;
        if(!is_euler)
            *X = g.normale();
        else{
            euler scheme(o,m,g);
            X = scheme.solution(EULER_PRECISION,true, 0.01);
        }
        mv -> maj( func(o,m,X,pas(compteur)) );
    }
    return compteur;
};

//Run une simulation et affiche le résultat sur stdout
void run( Model b, Option c, Generateur g, MeanVar *mv, std::function<double (int)> pas, bool is_euler, functions f ){
    std::clock_t start;
    double duration;
    start = std::clock();

    int compteur = MonteCarlo(c, b, g, f, pas, mv, is_euler);
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<< "S = " << b.S() << ", mean = "<< mv->mean() << ", var = " << mv -> var() <<", compteur = " << compteur << ", time = " << duration <<  std::endl;

}

int main(){

    double S = 100.0 ;
    double r = 0.02;
    double vol = 0.3;
    double nu = 0.7;
    //CHOIX DU MODELE

    BSmodel b(r, vol, S); 
    //CEVModel b(r,S,nu);

    bool is_euler = false; //mettre a true pour le model CEV

    //Paramètres Option
    double K = 100.0;
    double T = 5;
    Call c(K,T);

    /*Halton h(2); 
    Generateur g = LowDiscrepancyGenerator(h);*/

    Generateur g = RandomGenerator();

    //Pas de temps possible pour différence finie
   std::function<double (int)> pas_const =  [&] (int i){ return 0.001;};
   std::function<double (int)> pas_dcr = [&] (int i){ 
       return 1 / pow(i,0.25);
    };

    //fonction à tester (delta par malliavin, ou par différence finie..)
    //cf utils.hpp
    functions f  = Tangent_delta_bs;

    double precision = 0.001;
    int compteur;
    
    MeanVar *mv = new MeanVar;
    *mv = MeanVar(precision);
    //Lancer une simulation
    run(b,c,g,mv,pas_const, is_euler,f);

    //Où bien tracer une courbe complète

 /*   char const *output_name = "output.csv";
    std::ofstream fdm_out = std::ofstream(output_name);

    fdm_out << "S" << ";" << "mean" << ";" << "var"<< ";" << "nb_sim" << ";" << "time" << std::endl;

  /*  for(S = 10; S < 200 ; S+=2){
        std::clock_t start;
        double duration;
        start = std::clock();

        //Re-préciser le modèle ici
        BSmodel b(r,vol,S);
        //CEVModel b(r,S,nu);

        int compteur = MonteCarlo(c, b, g, f, pas_dcr, mv, is_euler);
        duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

        std::cout << "rate: " << b.rate(0.0,0.0) << ", vol: " << b.sigma(0.0,0.0) << ", strike: "<< c.Strike() << ", horizon: " << c.Horizon() << std::endl;
        std::cout<< "S = " << S << ", mean = "<< mv->mean() << ", var = " << mv -> var() <<", compteur = " << compteur << ", time = " << duration <<  std::endl;
        fdm_out << S << ";" << mv->mean() << ";" << mv->var() << ";" << compteur << ";" << duration << std::endl;
        *mv = MeanVar(precision);
    }
    fdm_out.close();*/
    delete mv;

    return 0;
}