#include "Halton.hpp"
#include <cmath>


Halton::Halton(long _dim, long index):
    LowDiscrepancySequence( [this] () { 
        std::vector<double> result;
        result.resize(dim, 0.0);
        for(long i = 0; i < dim; i++){
            result[i] = phi(p[i]);
        }
        incr();
        //std::cout << result[0] << " " << result[1] << " " << get_index() << std::endl;
        return result;
    } ,_dim),
    index(index)
{
    find_primes(_dim);
}

void Halton::incr(){
    index++;
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

double Halton::phi(long p){
    std::vector<long> a = decompose(index,p);
    double result = 0.0;
    for(unsigned int i = 0; i < a.size(); i++){
        result += a[i] / pow(p,i+1);
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



std::vector<long> Halton::decompose(long n, long p){
    std::vector<long> result;
    long r = 1;
    while(pow(p,r) <= n){
        r++;
    }
    result.resize(r,0);
    long a, power;
    for(long i = r ; i > 0 ; i--){
        power = static_cast<long>(pow(p,i-1));
        if(n / power > 0)
            a = n / power;
        else if(i > 1)
            a=0;
        else 
            a=n;
        result[i-1] = a;
        n -= a*power;
    }
    return result;
}