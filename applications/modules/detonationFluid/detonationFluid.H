/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | detonationFluid for OpenFOAM Foundation v14
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM and detonationFoam V2.0.
    It is distributed under the GNU General Public License version 3 or later.

Class
    Foam::solvers::detonationFluid

Description
    Density-based reacting multicomponent detonation solver module for foamRun.

    Gate-B5 implementation combines the OpenFOAM-14 shockFluid density-based
    framework with OpenFOAM-14 multicomponent chemistry/transport and restores
    OF14-native HLL, HLLC, HLLCP, AUSM+ and AUSM+up backends from the original
    detonationFoam formulation. Kurganov and Tadmor remain available and
    Kurganov remains the regression baseline.

    Supported solverType values in constant/solverTypeProperties:
      - Euler
      - NS_Sutherland

    NS_mixtureAverage from the original OpenFOAM-8 code is intentionally
    deferred until its custom transport model is migrated/qualified.
\*---------------------------------------------------------------------------*/

#ifndef detonationFluid_H
#define detonationFluid_H

#include "basicFluidSolver.H"
#include "fluidMulticomponentThermo.H"
#include "compressibleMomentumTransportModel.H"
#include "reactionModel.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "multivariateScheme.H"

namespace Foam
{
namespace solvers
{

class detonationFluid
:
    public basicFluidSolver
{
protected:

    // Thermophysical properties
    autoPtr<fluidMulticomponentThermo> thermoPtr_;
    fluidMulticomponentThermo& thermo_;
    volScalarField& p_;
    volScalarField rho_;

    // Kinematic/conservative fields
    volVectorField U_;
    surfaceScalarField phi_;
    volScalarField K_;

    // Composition
    PtrList<volScalarField>& Y_;
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields_;

    // Physics selection and diagnostics
    IOdictionary solverTypeProperties_;
    word solverType_;
    bool inviscid_;
    scalar shockPositionLimit_;
    scalar shockDetectionPressure_;
    Switch stopAtShockPosition_;
    scalar shockPosition_;

    volScalarField magGradrho_;
    volScalarField maxp_;
    volScalarField Qdot_;
    volScalarField gamma_;
    volScalarField Rgas_;

    // Native OF14 transport and reaction models
    autoPtr<compressible::momentumTransportModel> momentumTransport_;
    autoPtr<fluidMulticomponentThermophysicalTransportModel>
        thermophysicalTransport_;
    autoPtr<reactionModel> reaction_;

    // Central-upwind controls
    word fluxScheme_;

    // Cached face fields (same structure as OF14 shockFluid)
    tmp<surfaceScalarField> pos_;
    tmp<surfaceScalarField> neg_;
    tmp<surfaceScalarField> rhoPos_;
    tmp<surfaceScalarField> rhoNeg_;
    tmp<surfaceVectorField> rhoUPos_;
    tmp<surfaceVectorField> rhoUNeg_;
    tmp<surfaceVectorField> UPos_;
    tmp<surfaceVectorField> UNeg_;
    tmp<surfaceScalarField> pPos_;
    tmp<surfaceScalarField> pNeg_;
    tmp<surfaceScalarField> aPos_;
    tmp<surfaceScalarField> aNeg_;
    tmp<surfaceScalarField> aSf_;
    tmp<surfaceScalarField> aphivPos_;
    tmp<surfaceScalarField> aphivNeg_;

    // Unified conservative-face terms.  These reduce exactly to the native
    // central-upwind expressions for Kurganov/Tadmor/HLL and carry the HLLC
    // contact correction without introducing a parallel solver path.
    tmp<surfaceScalarField> pFlux_;
    tmp<surfaceScalarField> pEnergyFlux_;
    tmp<surfaceScalarField> pWork_;
    tmp<surfaceVectorField> Uf_;
    tmp<surfaceScalarField> amaxSf_;

    // HLLCP pressure-difference corrections, zero/invalid for all other
    // backends.  These are separate because the HLLCP mass correction is not
    // proportional to either reconstructed density state and its momentum
    // correction is generally not parallel to the face normal.
    tmp<surfaceScalarField> rhoPhiCorrection_;
    tmp<surfaceVectorField> rhoUPhiCorrection_;
    tmp<surfaceScalarField> rhoEPhiCorrection_;

    tmp<surfaceVectorField> devTau_;

    // Optional LTS reciprocal time-step field
    tmp<volScalarField> trDeltaT_;

private:

    void setRDeltaT(const surfaceScalarField& amaxSf);
    void correctAcousticCoNum(const surfaceScalarField& amaxSf);
    void correctDensity();
    void fluxPredictor();
    void clearTemporaryFields();
    void updateDiagnostics();

    template<class Type>
    inline tmp<SurfaceField<Type>> interpolate
    (
        const VolField<Type>& vf,
        const surfaceScalarField& dir,
        const word& reconFieldName = word::null
    )
    {
        tmp<SurfaceField<Type>> tsf
        (
            fvc::interpolate
            (
                vf,
                dir,
                "reconstruct("
              + (reconFieldName != word::null ? reconFieldName : vf.name())
              + ')'
            )
        );

        SurfaceField<Type>& sf = tsf.ref();
        sf.rename(vf.name() + '_' + dir.name());
        return tsf;
    }

public:

    // Public references for function objects and other framework components
    const fluidMulticomponentThermo& thermo;
    const volScalarField& p;
    const volScalarField& rho;
    const volVectorField& U;
    const surfaceScalarField& phi;
    const PtrList<volScalarField>& Y;

    TypeName("detonationFluid");

    explicit detonationFluid(fvMesh& mesh);
    detonationFluid(const detonationFluid&) = delete;
    virtual ~detonationFluid();

    virtual void preSolve();
    virtual void moveMesh();
    virtual void motionCorrector();
    virtual void prePredictor();
    virtual void momentumTransportPredictor();
    virtual void thermophysicalTransportPredictor();
    virtual void momentumPredictor();
    virtual void thermophysicalPredictor();
    virtual void pressureCorrector();
    virtual void momentumTransportCorrector();
    virtual void thermophysicalTransportCorrector();
    virtual void postSolve();

    void operator=(const detonationFluid&) = delete;
};

} // End namespace solvers
} // End namespace Foam

#endif

// ************************************************************************* //
