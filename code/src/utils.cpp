/**
 * Contient toutes les fonctions utilisées pour simuler les grecques
 * Malliavin, Différence finie, Processus tangent 
 */ 

#include "utils.hpp"

double mv_delta_bs(Option x, Model y, double G , double e){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
   // std::cout << " g= "<< G << " horizon = " << S   << std::endl;
    return exp(-y.rate(0.0,0.0)*x.Horizon())*x.payoff(S)*G*sqrt(x.Horizon())/(y.S()*y.sigma(0.0,0.0)*x.Horizon());
};

double mv_gamma_bs(Option x, Model y, double G , double e){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    return exp(-y.rate(0.0,0.0)*x.Horizon())*x.payoff(S)*(pow(G*sqrt(x.Horizon()),2)/(y.sigma(0.0,0.0)*x.Horizon())-G*sqrt(x.Horizon())-1/y.sigma(0.0,0.0))/(pow(y.S(),2)*y.sigma(0.0,0.0)*x.Horizon());
};

double mv_vega_bs(Option x, Model y, double G, double e ){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    return exp(-y.rate(0.0,0.0)*x.Horizon()) * x.payoff(S) * ( pow(G*sqrt(x.Horizon()),2) / (y.sigma(0.0,0.0)*x.Horizon()) -G*sqrt(x.Horizon()) - 1/y.sigma(0.0,0.0));
};


double Tangent_delta_bs(Option x, Model y, double G, double e ){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
   // std::cout << S << std::endl;
    if(S>=x.Strike()) {
        return exp(-y.rate(0.0,0.0)*x.Horizon())*S/y.S();
    } else {
        return 0;
    }
};

double Tangent_vega_bs(Option x, Model y, double G, double e ){
    double S = y.S()*exp((y.rate(0.0,0.0)-0.5*pow(y.sigma(0.0,0.0),2))*x.Horizon()+G*sqrt(x.Horizon())*y.sigma(0.0,0.0));
    if(S>=x.Strike()) {
        return exp(-y.rate(0.0,0.0)*x.Horizon()) * (G*sqrt(x.Horizon())-y.sigma(0.0,0.0)*x.Horizon())*S;
    }else 
       return 0;
};



double Tangent_delta_cev(Option x, Model y, double G, double e ){
    double S =G;
    if(S>=x.Strike()) {
        return exp(-y.rate(0.0,0.0)*x.Horizon())*S/y.S();
    } else {
        return 0;
    }
};

static double FD_delta_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0);
    double r = m.rate(0.0,0.0);


    double S1 = o.payoff( (S0+e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x*sqrt(o.Horizon()) ) );
    double S2 = o.payoff( (S0-e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x*sqrt(o.Horizon()) ) );
    return (S1 - S2)/(2*e);
}

static double FD_delta_cev(Option o, Model m, double x, double e){
    double S0 = m.S();

    double S1 = o.payoff( (S0+e) *  x /m.S() );
    double S2 = o.payoff( (S0-e) *  x /m.S() );
    return (S1 - S2)/(2*e);
} //FYI tu penses pas qu'il faut actualiser????????



static double FD_vega_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma1 = m.sigma(0.0,0.0) + e;
    double sigma2=m.sigma(0.0,0.0) -e;
    double r = m.rate(0.0,0.0);

    double S1 = o.payoff( (S0) *  exp( (r - 0.5*sigma1*sigma1) * o.Horizon() + (sigma1)*x*sqrt(o.Horizon()) ) );
    double S2 = o.payoff(  (S0) *  exp( (r - 0.5*sigma2*sigma2) * o.Horizon() + sigma2*x*sqrt(o.Horizon()) ) );
    return (S1 - S2)/(2*e);
}

static double FD_gamma_bs(Option o, Model m, double x, double e){
    double S0 = m.S();
    double sigma = m.sigma(0.0,0.0);
    double r = m.rate(0.0,0.0);
    //printf("\nFD_gamma_bs x = %f S0 = %f vol = %f r = %f \n",x,S0,sigma, r);fflush(stdout);


    double S1 = o.payoff( (S0+e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x*sqrt(o.Horizon()) ) );
    double S2 = o.payoff(  (S0-e) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x*sqrt(o.Horizon()) ) );
    double S3 = o.payoff( (S0) *  exp( (r - 0.5*sigma*sigma) * o.Horizon() + sigma*x*sqrt(o.Horizon()) ) );
   // std::cout << "\nFD_gamma_bs S1 / S2 / S3 :"<< S1 <<"/"<<S2 << "/" << S3 << std::endl;
   // std::cout << "\nresult = " << (S1 + S2 - 2*S3) / (e*e) << std::endl;

    return (S1 + S2 - 2*S3) / (e*e);
}

static double FD_gamma_cev(Option o, Model m, double x, double e){
    double S0 = m.S();
    //printf("\nFD_gamma_bs x = %f S0 = %f vol = %f r = %f \n",x,S0,sigma, r);fflush(stdout);
    
    
    double S1 = o.payoff( (S0+e) *  x /m.S() );
    double S2 = o.payoff(  (S0-e) * x /m.S() );
    double S3 = o.payoff(x);
    // std::cout << "\nFD_gamma_bs S1 / S2 / S3 :"<< S1 <<"/"<<S2 << "/" << S3 << std::endl;
    // std::cout << "\nresult = " << (S1 + S2 - 2*S3) / (e*e) << std::endl;
    
    return (S1 + S2 - 2*S3) / (e*e);
}




