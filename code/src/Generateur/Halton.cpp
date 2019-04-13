#include "Halton.hpp"
#include <cmath>


Halton::Halton(long _dim, long index): LowDiscrepancySequence(_dim), index(index){
    find_primes(_dim);
}

void Halton::find_primes(long _dim){
    std::vector<long> result;
    result.resize(_dim,0);
    long count = 0;
    long potential_prime = 2;
    while(count < _dim){
        if(is_prime(potential_prime)){
            result[count] = potential_prime;
            count++;
        }
        potential_prime++;
    }
    p = result;
}

std::vector<double> Halton::next(){
    std::vector<double> result;
    result.resize(dim, 0.0);

    for(long i = 0; i < dim; i++){
        result[i] = phi(p[i]);
    }
    index++;
    return result;
}

double Halton::phi(long p){
    std::vector<long> a = decompose(index,p);
    double result = 0.0;
    for(unsigned int i = 0; i < a.size(); i++){
       // printf( "   ai = %ld  phi pow = %f   ",a[i], pow(p,i+1));
        result += a[i] / pow(p,i+1);
       // printf("phi result %f ", result);
    }
    return result;
}

bool Halton::is_prime(long p){
    for(long i = 2; i <= p / 2; ++i){
      if(p % i == 0){
          return false;
      }
    }
    return true;
}


//A FINIR !!!!
std::vector<long> Halton::decompose(long n, long p){
   // printf(" Decompose %ld avec %ld : ",n,p);
    std::vector<long> result;
    long r = 1;
    while(pow(p,r) <= n){
        r++;
    }
    //printf(" pow = %f ", pow(p,r));
    result.resize(r,0);
    long a, power;
    for(long i = r ; i > 0 ; i--){
      //  printf(" p = %ld , i = %ld, ",p,i);
        power = static_cast<long>(pow(p,i-1));
      //  printf("power = %ld, n / power = %ld ", power,n / power);
        if(n / power > 0)
            a = n / power;
        else if(i > 1)
            a=0;
        else 
            a=n;
        result[i-1] = a;
        n -= a*power;
    }
  //  printf("result = ");
    for(int i = 0 ; i < result.size(); i++){
    //    printf(" %ld", result[i]);
    }
  //  printf("\n");
    return result;
}