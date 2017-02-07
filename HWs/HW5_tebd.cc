#include "itensor/all.h"

using namespace itensor;
using std::vector;

vector<ITensor>
makeH(SiteSet const& sites)
    {
    auto N = sites.N();
    auto H = vector<ITensor>(N+1);
    for(auto b : range1(N-1))
        {
        //Make S.S operator on sites b, b+1
        H.at(b) = sites.op("Sz",b)*sites.op("Sz",b+1)
                   + 0.5*sites.op("S+",b)*sites.op("S-",b+1)
                   + 0.5*sites.op("S-",b)*sites.op("S+",b+1);
        }
    return H;
    }

vector<ITensor>
makeGates(vector<ITensor> const& H,
          Real tau)
    {
    auto gates = H;
    for(auto& g : gates)
        {
        if(!g) continue;
        g = expHermitian(-tau*g);
        }
    return gates;
    }

void
doSVD(MPS & psi,
      int b,
      ITensor phi,
      Direction dir,
      Real cutoff)
    {
    auto U = psi.A(b);
    ITensor D,V;
    svd(phi,U,D,V,{"Cutoff",cutoff});
    if(dir == Fromleft)
        {
        //multiply D into V
        psi.setA(b,U);
        psi.setA(b+1,D*V);
        }
    else
        {
        //multiply D into U
        psi.setA(b,U*D);
        psi.setA(b+1,V);
        }
    }

ITensor
applyGate(ITensor phi,
          ITensor gate)
    {
    //TODO:
    //1. Apply gate to phi using * operator

    //2. Restore original prime level of phi's indices

    //3. Normalize phi by dividing by norm(phi)

    return phi;
    }


int 
main()
    {
    int N = 10;
    Real tau = 0.005;
    Real cutoff = 1E-12;
    int nsweep = 10000;

    auto sites = SpinHalf(N);

    auto state = InitState(sites);
    for(auto n : range1(N))
        {
        if(n%2==1) state.set(n,"Up");
        else       state.set(n,"Dn");
        }

    auto psi = MPS(state);

    auto H = makeH(sites);
    auto gates = makeGates(H,tau);

    for(auto sw : range1(nsweep))
        {
        println("Starting sweep ",sw);
        auto Energy = 0.0;
        for(auto b : range1(N-1))
            {
            auto phi = psi.A(b)*psi.A(b+1);

            //TODO: add measurement of energy <S.S>

            phi = applyGate(phi,gates.at(b));

            doSVD(psi,b,phi,Fromleft,cutoff);
            }
        printfln("  Half 1: Energy per site = %.10f",Energy/N);
        Energy = 0.0;
        for(auto b = N-1; b >= 1; --b)
            {
            auto phi = psi.A(b)*psi.A(b+1);

            //TODO: add measurement of energy <S.S>

            phi = applyGate(phi,gates.at(b));

            doSVD(psi,b,phi,Fromright,cutoff);

            if(b == N/2) 
                {
                println("Bond dimension at center = ",
                         commonIndex(psi.A(b),psi.A(b+1)).m());
                }
            }
        printfln("  Half 2: Energy per site = %.10f",Energy/N);
        }

    return 0;
    }
