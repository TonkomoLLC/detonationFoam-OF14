/*---------------------------------------------------------------------------*\
  Qualification-only transport probe for Candidate-1 published OF8 H/H2
  thermal diffusion.  The case is configured with uniform composition and a
  temperature gradient.  Only an H2-X thermal-diffusion polynomial is nonzero.
\*---------------------------------------------------------------------------*/
#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "fluidMulticomponentThermo.H"
#include "compressibleMomentumTransportModel.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "surfaceInterpolate.H"
#include "PstreamReduceOps.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    autoPtr<fluidMulticomponentThermo> thermoPtr
    (
        fluidMulticomponentThermo::New(mesh)
    );
    fluidMulticomponentThermo& thermo = thermoPtr();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        thermo.renameRho()
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensions::velocity
    );

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        linearInterpolate(rho*U) & mesh.Sf()
    );

    autoPtr<compressible::momentumTransportModel> momentumTransport
    (
        compressible::momentumTransportModel::New(rho, U, phi, thermo)
    );

    autoPtr<fluidMulticomponentThermophysicalTransportModel>
        thermophysicalTransport
        (
            fluidMulticomponentThermophysicalTransportModel::New
            (
                momentumTransport(),
                thermo
            )
        );

    const PtrList<volScalarField>& Y = thermo.Y();
    const volScalarField& YH = thermo.Y("H");
    const volScalarField& YH2 = thermo.Y("H2");

    const tmp<surfaceScalarField> tjH = thermophysicalTransport->j(YH);
    const tmp<surfaceScalarField> tjH2 = thermophysicalTransport->j(YH2);
    const surfaceScalarField& jH = tjH();
    const surfaceScalarField& jH2 = tjH2();

    surfaceScalarField sumJ
    (
        IOobject
        (
            "sumJQualification",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(jH.dimensions(), 0)
    );

    forAll(Y, i)
    {
        sumJ += thermophysicalTransport->j(Y[i]);
    }

    const surfaceScalarField YHf(linearInterpolate(YH));
    const surfaceScalarField YH2f(linearInterpolate(YH2));

    scalar maxAbsH = 0;
    scalar maxAbsH2 = 0;
    scalar maxAbsSum = 0;
    scalar maxAbsRelation = 0;
    scalar maxRelationScale = VSMALL;
    scalar maxAbsDirectH2Inference = 0;
    scalar maxDirectScale = VSMALL;
    label nFaces = 0;

    // On a one-dimensional slab the internal faces are the axial faces.
    // For the published OF8 assignment with uniform composition:
    //   raw_H  = q  (contains both H and H2 polynomials)
    //   raw_H2 = 0
    //   sumRaw = q
    // hence
    //   j_H  = -(1-Y_H) q
    //   j_H2 = Y_H2 q
    // and therefore (1-Y_H) j_H2 + Y_H2 j_H = 0.
    // The inferred direct H2 raw contribution is
    //   raw_H2 = Y_H2*q - j_H2,
    // with q=-j_H/(1-Y_H), and should be zero.
    for (label facei=0; facei<mesh.nInternalFaces(); ++facei)
    {
        const scalar yh = YHf[facei];
        const scalar yh2 = YH2f[facei];
        const scalar h = jH[facei];
        const scalar h2 = jH2[facei];
        const scalar sj = sumJ[facei];

        const scalar relation = (1.0-yh)*h2 + yh2*h;
        const scalar relationScale =
            max(max(mag((1.0-yh)*h2), mag(yh2*h)), VSMALL);

        const scalar qFromH = -h/max(1.0-yh, VSMALL);
        const scalar inferredRawH2 = yh2*qFromH - h2;
        const scalar directScale = max(max(mag(yh2*qFromH), mag(h2)), VSMALL);

        maxAbsH = max(maxAbsH, mag(h));
        maxAbsH2 = max(maxAbsH2, mag(h2));
        maxAbsSum = max(maxAbsSum, mag(sj));
        maxAbsRelation = max(maxAbsRelation, mag(relation));
        maxRelationScale = max(maxRelationScale, relationScale);
        maxAbsDirectH2Inference =
            max(maxAbsDirectH2Inference, mag(inferredRawH2));
        maxDirectScale = max(maxDirectScale, directScale);
        ++nFaces;
    }

    reduce(maxAbsH, maxOp<scalar>());
    reduce(maxAbsH2, maxOp<scalar>());
    reduce(maxAbsSum, maxOp<scalar>());
    reduce(maxAbsRelation, maxOp<scalar>());
    reduce(maxRelationScale, maxOp<scalar>());
    reduce(maxAbsDirectH2Inference, maxOp<scalar>());
    reduce(maxDirectScale, maxOp<scalar>());
    reduce(nFaces, sumOp<label>());

    const scalar relationRel = maxAbsRelation/maxRelationScale;
    const scalar inferredDirectH2Rel = maxAbsDirectH2Inference/maxDirectScale;
    const scalar fluxScale = max(max(maxAbsH,maxAbsH2), VSMALL);
    const scalar sumRel = maxAbsSum/fluxScale;

    Info<< "SORET_METRIC internalFaces=" << nFaces
        << " maxAbsJH=" << maxAbsH
        << " maxAbsJH2=" << maxAbsH2
        << " maxAbsSumJ=" << maxAbsSum
        << " sumRel=" << sumRel
        << " publishedRelationRel=" << relationRel
        << " inferredDirectH2Rel=" << inferredDirectH2Rel
        << endl;
    Info<< "SORET_END" << endl;

    return 0;
}
