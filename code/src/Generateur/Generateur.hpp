#ifndef GENERATEUR_H
#define GENERATEUR_H

#include <functional>

class Generateur{

    public:
        Generateur(std::function<double ()> _normale): normale(_normale){};

        std::function<double ()> normale;
};

#endif