#include "detonationFluid.H"

void Foam::solvers::detonationFluid::moveMesh()
{
    if (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    {
        mesh_.move();
    }
}


void Foam::solvers::detonationFluid::motionCorrector()
{
    if (pimple.firstIter() || pimple.moveMeshOuterCorrectors())
    {
        if (mesh.changing())
        {
            meshCourantNo();
        }
    }
}

// ************************************************************************* //
