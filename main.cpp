#include <vector>
#include <iostream>
#include "bulk/bulk.hpp"
#include "../examples/set_backend.hpp"
#include "includes/fieldArithmetic.h"

//I think we don't need fwd or n like in the regular algorithm, since fwd is encoded in the weight table.
void butterflyStage(std::vector<unsigned long> & x, nntlib::finField & field, unsigned long n, unsigned long k,
                    std::vector<unsigned long> & w, unsigned long start){
    for (long r = 0; r < n/k; r++){
        for (long j = 0; j < k/2; j++){
            //Like discussed, just the same using field arithmetic.
            unsigned long tau = field.mult(w[j + start],x[r*k+j+k/2]);
            x[r*k+j+k/2] = field.add(x[r*k+j], field.addInv(tau));
            x[r*k+j] = field.add(x[r*k+j], tau);
        }
    }
}

void unnt(std::vector<unsigned long> & x, nntlib::finField & field, unsigned long n, std::vector<unsigned long> & w){
    unsigned long start = 0;
    for (long k = 2; k <= n; k *= 2){
        butterflyStage(x, field, n, k, w, start); //Need the start since want to pass where the weight table should start.
        start += k/2;
    }
}

void permute(std::vector<unsigned long> & x, unsigned long n, std::vector<unsigned long> & sigm){
    for (unsigned long j = 0; j < n; j++){
        unsigned long sigmaj = sigm[j];
        if (j < sigmaj) {
            unsigned long temp = x[j];
            x[j] = x[sigmaj];
            x[sigmaj] = temp;
        }
    }
}

void bitrev_init(unsigned long n, std::vector<unsigned long> & rho){
    rho[0] = 0;
    for (unsigned long k = 2; k <= n; k *= 2){
        for (unsigned long j = 0; j < k/2; j++){
            rho[j] *= 2;
            rho[j+k/2] = rho[j]+1;
        }
    }
}

void bspredistr(bulk::world & world, std::vector<unsigned long> & x, unsigned long c0, unsigned long c, bool rev, std::vector<unsigned long> & rho_p){
    long s = world.rank();
    long np = x.size();

    long j0, j2;
    if (rev){
        j0 = rho_p[s] % c0;
        j2 = rho_p[s] / c0;
    } else {
        j0 = s % c0;
        j2 = s / c0;
    }
    long ratio = c/c0;
    long size = (np >= ratio ? np/ratio : 1);
    long npackets = np/size;

    //Using a queue instead of regular put, since the exact memory reference isnt the same. Still equivalent.
    auto q = bulk::queue<unsigned long, unsigned long>(world);

    for (long j = 0; j < npackets; j++){
        long jglob = j2*c0*np + j*c0 + j0;
        long destproc = (jglob/(c*np))*c + jglob%c;
        long destindex = (jglob %(c*np)) / c;
        for (long r = 0; r < size; r++){
            q(destproc).send(destindex+r, x[j+r*ratio]);
        }
    }
    world.sync();
    for (auto [ind, val] : q){
        x[ind]=val;
    }
}

void bspnnt(bulk::world & world, std::vector<unsigned long> & x, std::vector<unsigned long> & w,
            std::vector<unsigned long> & rho_np, std::vector<unsigned long> & rho_p, nntlib::finField & field, bool fwd){
    unsigned long p = world.active_processors();
    unsigned long np = x.size();
    unsigned long c = 1;
    bool rev = true;

    permute(x, np, rho_np);
    unnt(x, field, np, w);

    unsigned long k = 2*np;
    unsigned long start = np-1;
    while(c < p){
        unsigned long c0 = c;
        c = (np*c <= p ? np*c : p);
        bspredistr(world, x,c0, c ,rev,rho_p);
        rev = false;

        while(k <= np*c){
            butterflyStage(x, field, np, k/c, w, start);
            start += k/(2*c);
            k*= 2;
        }
    }
    /*normalise?*/
    if (!fwd){
        for (long i = 0 ; i < x.size(); i++)
            x[i] = field.divbypow2(x[i], np * p);
    }
}

void bspnnt_init(bulk::world & world, unsigned long n, std::vector<unsigned long> & w, nntlib::finField & field,
                 std::vector<unsigned long> & rho_np, std::vector<unsigned long> & rho_p, bool fwd){
    long s = world.rank();
    long p = world.active_processors();
    long np = n/p;

    //bitrev perms
    bitrev_init(np, rho_np);
    bitrev_init(p, rho_p);

    //Field

    // omega tables
    long c = 1;
    //long k = 2;
    unsigned long kxpo = 1; // 2^(k) == 1<<k.
    // so ((unsigned long)1)<<kxpo == k in book.
    unsigned long w0 = (fwd ? 182099332568824125 : 68630377364883);
    long start = 0;
    while (c <= p){
        unsigned long j0 = s%c;
        while ((((unsigned long)1)<<kxpo) <= np*c){
            unsigned long basew = field.pow(w0, (((unsigned long)1)<<(57-kxpo)));
            // In fft, we use -2pi / k, which represents a k-th root. I think this is equivalent to that
            unsigned long cfactor = field.pow(basew, c);
            unsigned long wcalc = field.pow(basew, j0);
            // in fft each forloop, we calculate the exponent as e^ (i*(j0+j*c)*theta), so we add c each time except for j=0;
            // Here we want w^((2^57/k)*(j0+j*c)). = wcalc*(cfactor)^j = (w^(2^(2^57/k)))^j0 * (w^(2^57*c))^j.
            if (0 < (((unsigned long)1)<<kxpo)/(2*c)) {
                w[start] = wcalc;
                start++;
            }
            for (long j = 1; j < (((unsigned long)1)<<kxpo)/(2*c); j++, start++){
                //With fft, we don't want to use previous weights for new ones due to accuracy, but in ntt we can.
                wcalc = field.mult(wcalc, cfactor);
                w[start] = wcalc;
            }
            kxpo++;
        }
        if (c < p)
            // if "c small enough"  then c = np*c.
            c = (np*c <= p ? np*c : p);
        else
            c = 2*p;
    }
}

void nntcall_p(bulk::world & world, std::vector<unsigned long> & localx, bool fwd){
    long p = world.active_processors();
    long np = localx.size();
    long n = np * p;
    long s = world.rank();

    std::vector<unsigned long> rho_np(np);
    std::vector<unsigned long> rho_p(p);
    std::vector<unsigned long> w(n);
    auto field = nntlib::finField(4179340454199820289);

    // DEBUGGING values
    int debugproc = 0;
    if (s == debugproc){
        std::cout << "s: " << s << " - " << "[";
        int counter=0;
        for (auto elem : localx){
            std::cout << elem << ", ";
            counter++;
            if (counter==500)
                break;
        }
        std::cout << "]" << std::endl;
    }
    // */// END DEBUGGING goldbach


    bspnnt_init(world, n, w, field, rho_np, rho_p, fwd);
    bspnnt(world, localx, w, rho_np, rho_p, field, fwd);

    // DEBUGGING values
    int counter=0;
    if (s == debugproc){
        std::cout << "s: " << s << " - " << "[";
        for (auto elem : localx){
            std::cout << elem << ", ";
            counter++;
            if (counter==500)
                break;
        }
        std::cout << "]" << std::endl;
    }
    // */// END DEBUGGING goldbach
    world.sync();
}

void nntcall(std::vector<unsigned long> & x, bool fwd){
    environment env;

    // I'm really, really not sure if this is the correct way to get the variables into the function, as this will not
    // distribute x.
    env.spawn(env.available_processors(), [x, fwd](bulk::world& world) {
        long p = world.active_processors();
        long n = x.size();
        long np = n/p;
        long s = world.rank();

        std::vector<unsigned long> localx(np);
        for (long i = 0; i < x.size(); i++) {
            if (i%p == s)
                localx[i/p] = x[i];
        }

        nntcall_p(world, localx, fwd);
    });
}

void nntcall_funct(std::function<unsigned long(unsigned long)> & f, unsigned long n,bool fwd){
    environment env;

    // I'm really, really not sure if this is the correct way to get the variables into the function, as this will not
    // distribute x.
    env.spawn(env.available_processors(), [f, n, fwd](bulk::world& world) {
        long p = world.active_processors();
        long np = n/p;
        long s = world.rank();

        std::vector<unsigned long> localx(np);
        for (long i = s; i < n; i+= p) {
            localx[i/p] = f(i) % 4179340454199820289;
        }

        nntcall_p(world, localx, fwd);
    });
}


int main() {

    int size;
    std::cin>>size;

    auto start = std::chrono::system_clock::now();

    /*
    std::vector<unsigned long> dummyX(16777216);
    for (int i = 0; i < 16777216; i++)
        dummyX[i] = 1;

    nntcall(dummyX, true);

    //*///

    auto f = std::function<unsigned long(unsigned long)> ([](unsigned long indexOfX){
        return indexOfX % 2;
    });
    nntcall_funct(f, 1<<size,true);
    //*///

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    //timer functionality from:
    // https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
    return 0;

}