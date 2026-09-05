/*---------------------------------------------------------------------------*\
  detonationFluid - OpenFOAM Foundation v14 foamRun solver module
\*---------------------------------------------------------------------------*/

#include "detonationFluid.H"
#include "fvMeshStitcher.H"
#include "localEulerDdtScheme.H"
#include "fvcGrad.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(detonationFluid, 0);
    addToRunTimeSelectionTable(solver, detonationFluid, fvMesh);
}
}

void Foam::solvers::detonationFluid::correctAcousticCoNum
(
    const surfaceScalarField& amaxSf
)
{
    const scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

    CoNum_ =
        0.5*gMax(sumAmaxSf/mesh.V().primitiveField())*runTime.deltaTValue();

    const scalar meanCoNum =
        0.5
       *(gSum(sumAmaxSf)/gSum(mesh.V().primitiveField()))
       *runTime.deltaTValue();

    Info<< "Acoustic Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}


void Foam::solvers::detonationFluid::clearTemporaryFields()
{
    rhoPos_.clear();
    rhoNeg_.clear();
    rhoUPos_.clear();
    rhoUNeg_.clear();
    UPos_.clear();
    UNeg_.clear();
    pPos_.clear();
    pNeg_.clear();
    aPos_.clear();
    aNeg_.clear();
    aSf_.clear();
    aphivPos_.clear();
    aphivNeg_.clear();
    pFlux_.clear();
    pEnergyFlux_.clear();
    pWork_.clear();
    Uf_.clear();
    amaxSf_.clear();
    rhoPhiCorrection_.clear();
    rhoUPhiCorrection_.clear();
    rhoEPhiCorrection_.clear();
    devTau_.clear();
}


Foam::solvers::detonationFluid::detonationFluid(fvMesh& mesh)
:
    basicFluidSolver(mesh),

    thermoPtr_(fluidMulticomponentThermo::New(mesh)),
    thermo_(thermoPtr_()),
    p_(thermo_.p()),

    rho_
    (
        IOobject
        (
            "rho",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo_.renameRho()
    ),

    U_
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensions::velocity
    ),

    phi_
    (
        IOobject
        (
            "phi",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho_*U_) & mesh.Sf()
    ),

    K_("K", 0.5*magSqr(U_)),

    Y_(thermo_.Y()),

    solverTypeProperties_
    (
        IOobject
        (
            "solverTypeProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    solverType_
    (
        solverTypeProperties_.lookupOrDefault<word>
        (
            "solverType",
            "NS_Sutherland"
        )
    ),

    inviscid_(solverType_ == "Euler"),

    shockPositionLimit_
    (
        solverTypeProperties_.lookupOrDefault<scalar>
        (
            "SW_position_limit",
            great
        )
    ),

    shockDetectionPressure_
    (
        solverTypeProperties_.lookupOrDefault<scalar>
        (
            "shockDetectionPressure",
            101325 + 100
        )
    ),

    stopAtShockPosition_
    (
        solverTypeProperties_.lookupOrDefault<Switch>
        (
            "stopAtShockPosition",
            true
        )
    ),

    shockPosition_(0),

    magGradrho_
    (
        IOobject
        (
            "magGradrho",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(fvc::grad(rho_))
    ),

    maxp_
    (
        IOobject
        (
            "maxp",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("maxp", dimPressure, 0)
    ),

    Qdot_
    (
        IOobject
        (
            "Qdot",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Qdot", dimEnergy/dimVolume/dimTime, 0)
    ),

    gamma_
    (
        IOobject
        (
            "gama",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo_.gamma()
    ),

    Rgas_
    (
        IOobject
        (
            "R_gas",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo_.Cp() - thermo_.Cv()
    ),

    momentumTransport_
    (
        compressible::momentumTransportModel::New
        (
            rho_,
            U_,
            phi_,
            thermo_
        )
    ),

    thermophysicalTransport_
    (
        fluidMulticomponentThermophysicalTransportModel::New
        (
            momentumTransport_(),
            thermo_
        )
    ),

    reaction_(reactionModel::New(thermo_, momentumTransport_())),

    fluxScheme_
    (
        mesh.schemes().lookupOrDefault<word>("fluxScheme", "Kurganov")
    ),

    thermo(thermo_),
    p(p_),
    rho(rho_),
    U(U_),
    phi(phi_),
    Y(Y_)
{
    thermo_.validate(type(), "e");
    momentumTransport_->validate();

    if (solverType_ != "Euler" && solverType_ != "NS_Sutherland")
    {
        FatalErrorInFunction
            << "Unsupported solverType '" << solverType_ << "'." << nl
            << "Gate A supports Euler and NS_Sutherland." << nl
            << "NS_mixtureAverage is intentionally deferred until the custom "
            << "mixture-average transport model is migrated and qualified."
            << exit(FatalError);
    }

    if
    (
        fluxScheme_ != "Kurganov"
     && fluxScheme_ != "Tadmor"
     && fluxScheme_ != "HLL"
     && fluxScheme_ != "HLLC"
     && fluxScheme_ != "HLLCP"
     && fluxScheme_ != "AUSM+"
     && fluxScheme_ != "AUSM+up"
    )
    {
        FatalErrorInFunction
            << "Unsupported Gate-B5 fluxScheme '" << fluxScheme_ << "'." << nl
            << "Supported schemes are Kurganov, Tadmor, HLL, HLLC, HLLCP, "
            << "AUSM+ and AUSM+up."
            << exit(FatalError);
    }

    Info<< "detonationFluid solverType = " << solverType_ << nl
        << "detonationFluid fluxScheme = " << fluxScheme_ << nl
        << "native OF14 chemistry/transport enabled" << nl
        << "legacy external chemistry and custom 2-D mesh libraries are not linked" << endl;

    forAll(Y_, i)
    {
        fields_.add(Y_[i]);
        Y_[i].writeOpt() = IOobject::AUTO_WRITE;
    }
    fields_.add(thermo_.he());

    mesh.schemes().setFluxRequired(U_.name());

    fluxPredictor();

    if (transient())
    {
        correctAcousticCoNum(amaxSf_());
    }
    else if (LTS)
    {
        Info<< "Using LTS" << endl;
        trDeltaT_ = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    fv::localEulerDdt::rDeltaTName,
                    runTime.name(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless/dimTime, 1),
                extrapolatedCalculatedFvPatchScalarField::typeName
            )
        );
    }
}


Foam::solvers::detonationFluid::~detonationFluid()
{}


void Foam::solvers::detonationFluid::preSolve()
{
    updateDiagnostics();

    {
        if (transient())
        {
            correctAcousticCoNum(amaxSf_());
        }
        else if (LTS)
        {
            setRDeltaT(amaxSf_());
        }
    }

    fvModels().preUpdateMesh();

    // Face-reconstruction temporaries are topology/addressing dependent.
    // fvMeshDistributor redistribution is not necessarily advertised through
    // topoChanging() before mesh_.update(); distributing() explicitly covers
    // that path. Clear them before redistribution so fluxPredictor() rebuilds
    // them against the new face/processor addressing.
    if (mesh.topoChanging() || mesh.distributing() || mesh.stitcher().stitches())
    {
        pos_.clear();
        neg_.clear();
        clearTemporaryFields();
    }

    // This is the OF14-native hook for topology change and mesh distribution.
    mesh_.update();
}


void Foam::solvers::detonationFluid::prePredictor()
{
    fluxPredictor();
    correctDensity();
}


void Foam::solvers::detonationFluid::momentumTransportPredictor()
{
    if (!inviscid_)
    {
        momentumTransport_->predict();
    }
}


void Foam::solvers::detonationFluid::thermophysicalTransportPredictor()
{
    if (!inviscid_)
    {
        thermophysicalTransport_->predict();
    }
}


void Foam::solvers::detonationFluid::momentumTransportCorrector()
{
    if (!inviscid_)
    {
        momentumTransport_->correct();
    }
}


void Foam::solvers::detonationFluid::thermophysicalTransportCorrector()
{
    if (!inviscid_)
    {
        thermophysicalTransport_->correct();
    }
}


void Foam::solvers::detonationFluid::postSolve()
{
    updateDiagnostics();

    Info<< "min/max(p) = " << min(p_).value() << ", " << max(p_).value()
        << nl
        << "min/max(T) = " << min(thermo_.T()).value() << ", "
        << max(thermo_.T()).value() << endl;

    scalar localShockPosition = -great;
    forAll(p_, celli)
    {
        if (p_[celli] > shockDetectionPressure_)
        {
            localShockPosition = max(localShockPosition, mesh.C()[celli].x());
        }
    }
    reduce(localShockPosition, maxOp<scalar>());

    if (localShockPosition > -0.5*great)
    {
        shockPosition_ = max(shockPosition_, localShockPosition);
    }

    Info<< "Leading shock location = " << shockPosition_ << " m" << endl;

    if (stopAtShockPosition_ && shockPosition_ >= shockPositionLimit_)
    {
        Info<< "Shock-position limit reached (" << shockPositionLimit_
            << " m); ending after current write." << endl;
        // solver::runTime is a const Time& in the OF14 module API.
        // Time controls are intentionally mutable at the application level;
        // use the same const_cast pattern as OF14 time-control functionObjects.
        const_cast<Time&>(runTime).setEndTime(runTime.value());
    }
}


void Foam::solvers::detonationFluid::updateDiagnostics()
{
    magGradrho_ = mag(fvc::grad(rho_));

    forAll(maxp_, celli)
    {
        if (p_[celli] > maxp_[celli])
        {
            maxp_[celli] = p_[celli];
        }
    }

    gamma_ = thermo_.gamma();
    Rgas_ = thermo_.Cp() - thermo_.Cv();
}

// ************************************************************************* //
