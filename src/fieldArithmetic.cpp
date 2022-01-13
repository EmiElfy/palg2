//
// Created by emi on 20/12/2021.
//
#include "../includes/fieldArithmetic.h"
#include <stdexcept>

namespace nntlib{
    finField::finField(unsigned long n) {
        if (n>>63 == 0)
            n_ = n;
        //else throw std::invalid_argument("input number is too large");
    }

    unsigned long finField::addInv(unsigned long a){
        if (a<=n_)
            return (n_ - a) % n_;
        //else throw std::invalid_argument("input number is too large");
    }

    unsigned long finField::add(unsigned long a, unsigned long b) {
        if (a < n_ && b < n_){
            return (a + b) % n_;
        }
        //else throw std::invalid_argument("one of the inputs is too large");
    }

    unsigned long finField::mult(unsigned long a, unsigned long b) {
        if (a < n_ && b < n_){
            __int128 anew = a;
            __int128 bnew = b;
            __int128 nnew = n_;
            return (anew * bnew) % nnew;
        }
        //else throw std::invalid_argument("one of the inputs is too large");
    }

    unsigned long finField::pow(unsigned long a, unsigned long b) {
        if (a < n_){
            unsigned long result = 1;
            while (true){
                if (b & 1)
                    result = mult(result, a);
                b >>= 1;
                if (!b) break;
                a = mult(a, a);
            }
            return result;
        }
        //else throw std::invalid_argument("one of the inputs is too large");
    }

    unsigned long finField::divbypow2(unsigned long a, unsigned long b){
        if (a < n_ && b < n_) {
            return addInv(mult(a, (n_ - 1) / b));
        }
    }
}