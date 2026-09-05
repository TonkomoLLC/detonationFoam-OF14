#include "detonationFluid.H"
#include "fvmDdt.H"
#include "fvcDiv.H"

void Foam::solvers::detonationFluid::correctDensity()
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho_) + fvc::div(phi_)
      ==
        fvModels().source(rho_)
    );

    fvConstraints().constrain(rhoEqn);
    rhoEqn.solve();
    fvConstraints().constrain(rho_);
}

// ************************************************************************* //
