#include "detonationFluid.H"
#include "fvcSmooth.H"

void Foam::solvers::detonationFluid::setRDeltaT
(
    const surfaceScalarField& amaxSf
)
{
    volScalarField& rDeltaT = trDeltaT_.ref();
    const dictionary& pimpleDict = pimple.dict();

    const scalar maxCo
    (
        pimpleDict.lookupOrDefault<scalar>("maxCo", 0.8)
    );
    const scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    rDeltaT.internalFieldRef() =
        fvc::surfaceSum(amaxSf)/((2*maxCo)*mesh.V());

    scalar minRDeltaT = gMin(rDeltaT.primitiveField());

    if (pimpleDict.found("maxDeltaT") || minRDeltaT < rootVSmall)
    {
        const scalar clipRDeltaT = 1/pimpleDict.lookup<scalar>("maxDeltaT");
        rDeltaT.max(clipRDeltaT);
        minRDeltaT = max(minRDeltaT, clipRDeltaT);
    }

    if (pimpleDict.found("minDeltaT"))
    {
        const scalar clipRDeltaT = 1/pimpleDict.lookup<scalar>("minDeltaT");
        rDeltaT.min(clipRDeltaT);
        minRDeltaT = min(minRDeltaT, clipRDeltaT);
    }

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField()) << ", " << 1/minRDeltaT << endl;

    rDeltaT.correctBoundaryConditions();
    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField()) << ", "
        << gMax(1/rDeltaT.primitiveField()) << endl;
}

// ************************************************************************* //
