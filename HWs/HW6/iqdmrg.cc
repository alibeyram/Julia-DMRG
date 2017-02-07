#include "itensor/all.h"

using namespace itensor;

int 
main(int argc, char* argv[])
    {
    int N = 100;

    //
    // Initialize the site degrees of freedom.
    //
    //auto sites = SpinHalf(N); //make a chain of N spin 1/2's
    auto sites = SpinOne(N); //make a chain of N spin 1's

    //
    // Use the AutoMPO feature to create the 
    // next-neighbor Heisenberg model.
    //
    // Here we convert the AutoMPO information
    // into an IQMPO, a matrix-product operator
    // which automatically tracks quantum
    // number information.
    //
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = IQMPO(ampo);

    //
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    // This choice implicitly sets the global Sz quantum number
    // of the wavefunction to zero. Since it is an IQMPS
    // it will remain in this quantum number sector.
    //
    //
    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i) 
        {
        if(i%2 == 1)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
        }

    state.set(2, "Up");  // flip one spin to go to exited state since itensor would preserve the total S_z (quantum number)
    // state.set(98, "Up");  // flip second spin somehow related to Haldane gap 

    auto psi = IQMPS(state);

    //
    // overlap(psi,H,psi) = <psi|H|psi>
    //
    printfln("Initial energy = %.5f", overlap(psi,H,psi) );

    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    // Here less than 5 cutoff values are provided, for example,
    // so all remaining sweeps will use the last one given (= 1E-10).
    //
    auto sweeps = Sweeps(30);
    sweeps.maxm() = 20,50,100,150,200,400;   // I changed them form 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;  
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    println(sweeps);

    //
    // Begin the DMRG calculation
    //
    auto energy = dmrg(psi,H,sweeps,"Quiet");

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing overlap = %.10f", overlap(psi,H,psi) );

    println("\nTotal QN of Ground State = ",totalQN(psi));

    return 0;
    }
