#include "detonationFluid.H"

void Foam::solvers::detonationFluid::pressureCorrector()
{
    const volScalarField& psi = thermo_.psi();
    p_.internalFieldRef() = rho_/psi;
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() == psi.boundaryField()*p_.boundaryField();
}

// ************************************************************************* //
