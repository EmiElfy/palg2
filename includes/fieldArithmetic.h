//
// Created by emi on 20/12/2021.
//

#ifndef TEMP_FIELDARITHMETIC_H
#define TEMP_FIELDARITHMETIC_H
namespace nntlib {
    class finField{
    public:
        finField(unsigned long n);
        unsigned long add(unsigned long a, unsigned long b);
        unsigned long mult(unsigned long a, unsigned long b);
        unsigned long addInv(unsigned long a);
        unsigned long pow(unsigned long a, unsigned long b);
        unsigned long divbypow2(unsigned long a, unsigned long b);
            private:
        unsigned long n_;
    };
}
#endif //TEMP_FIELDARITHMETIC_H
