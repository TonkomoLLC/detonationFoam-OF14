#include "detonationFluid.H"
#include "fvmDdt.H"
#include "fvcDiv.H"

void Foam::solvers::detonationFluid::momentumPredictor()
{
    surfaceVectorField phiUp
    (
        (aphivPos_()*rhoUPos_() + aphivNeg_()*rhoUNeg_())
      + pFlux_()*mesh.Sf()
    );

    if (rhoUPhiCorrection_.valid())
    {
        phiUp += rhoUPhiCorrection_();
    }

    tmp<fvVectorMatrix> divDevTau;
    if (!inviscid_)
    {
        divDevTau = momentumTransport_->divDevTau(U_);
    }

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho_, U_) + fvc::div(phiUp)
      ==
        fvModels().source(rho_, U_)
    );

    if (!inviscid_)
    {
        UEqn += divDevTau();
    }

    UEqn.relax();
    fvConstraints().constrain(UEqn);
    solve(UEqn);
    fvConstraints().constrain(U_);

    K_ = 0.5*magSqr(U_);

    if (!inviscid_)
    {
        devTau_ = divDevTau->flux();
    }
}

// ************************************************************************* //
