#include "legacyMixtureAverageFourier.H"
#include "Function2Evaluate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "surfaceInterpolate.H"

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

template<class BasicThermophysicalTransportModel>
const Function2<scalar>&
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::Dij
(
    const label i,
    const label j
) const
{
    return i < j ? DFuncs_[i][j] : DFuncs_[j][i];
}


template<class BasicThermophysicalTransportModel>
void legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::clearCache() const
{
    Dm_.clear();
    sumRawFlux_.clear();
}


template<class BasicThermophysicalTransportModel>
void legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::updateDm() const
{
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();
    const volScalarField Wm(this->thermo().W());

    Dm_.setSize(Y.size());

    forAll(Y, i)
    {
        volScalarField sumXbyD
        (
            volScalarField::New
            (
                "legacySumXbyD_" + Y[i].name(),
                T.mesh(),
                dimensionedScalar
                (
                    dimless/dimKinematicViscosity/Wm.dimensions(),
                    0
                )
            )
        );

        forAll(Y, j)
        {
            if (j != i)
            {
                sumXbyD +=
                    Y[j]
                   /(
                        this->thermo().Wi(j)
                       *evaluate
                        (
                            Dij(i, j),
                            dimKinematicViscosity,
                            p,
                            T
                        )
                    );
            }
        }

        Dm_.set
        (
            i,
            volScalarField::New
            (
                "legacyDm_" + Y[i].name(),
                ((scalar(1) - Y[i])/Wm)
               /max
                (
                    sumXbyD,
                    dimensionedScalar(sumXbyD.dimensions(), small)
                )
            ).ptr()
        );
    }
}


template<class BasicThermophysicalTransportModel>
const PtrList<volScalarField>&
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::Dm() const
{
    if (!Dm_.size())
    {
        updateDm();
    }
    return Dm_;
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::rawFlux
(
    const label i
) const
{
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField Wm(this->thermo().W());
    const tmp<volScalarField> tDi(this->momentumTransport().rho()*Dm()[i]);
    const volScalarField& Di = tDi();

    const tmp<volVectorField> tRaw
    (
        Di*fvc::grad(Y[i])
      + Di*Y[i]/Wm*fvc::grad(Wm)
    );

    return surfaceScalarField::New
    (
        "legacyRawFlux_" + Y[i].name(),
        (
            fvc::interpolate(this->alpha()*tRaw())
          & Y[i].mesh().Sf()
        )/Y[i].mesh().magSf()
    );
}


template<class BasicThermophysicalTransportModel>
void legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::updateSumRawFlux() const
{
    const PtrList<volScalarField>& Y = this->thermo().Y();
    sumRawFlux_ = surfaceScalarField::New
    (
        "legacySumRawFlux",
        Y[0].mesh(),
        dimensionedScalar(dimMass/dimArea/dimTime, 0)
    );
    forAll(Y, i)
    {
        sumRawFlux_.ref() += rawFlux(i);
    }
}


template<class BasicThermophysicalTransportModel>
const surfaceScalarField&
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::sumRawFlux() const
{
    if (!sumRawFlux_.valid())
    {
        updateSumRawFlux();
    }
    return sumRawFlux_();
}


template<class BasicThermophysicalTransportModel>
tmp<volScalarField>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::legacyKappa() const
{
    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();
    const volScalarField Wm(this->thermo().W());

    const tmp<volScalarField> tk0(this->thermo().kappai(0, p, T));
    const dimensionSet kappaDims(tk0().dimensions());

    volScalarField kappaArithmetic
    (
        volScalarField::New
        (
            "legacyKappaArithmetic",
            T.mesh(),
            dimensionedScalar(kappaDims, 0)
        )
    );

    volScalarField invKappaHarmonic
    (
        volScalarField::New
        (
            "legacyInvKappaHarmonic",
            T.mesh(),
            dimensionedScalar(dimless/kappaDims, 0)
        )
    );

    forAll(Y, i)
    {
        const volScalarField Xi(Y[i]*Wm/this->thermo().Wi(i));
        const tmp<volScalarField> tKappai(this->thermo().kappai(i, p, T));
        const volScalarField& kappai = tKappai();

        kappaArithmetic += Xi*kappai;
        invKappaHarmonic +=
            Xi/max(kappai, dimensionedScalar(kappaDims, small));
    }

    return volScalarField::New
    (
        "legacyKappa",
        scalar(0.5)
       *(
            kappaArithmetic
          + scalar(1)
           /max
            (
                invKappaHarmonic,
                dimensionedScalar(dimless/kappaDims, small)
            )
        )
    );
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::legacyKappa
(
    const label patchi
) const
{
    const tmp<volScalarField> tkappa(this->legacyKappa());
    return tmp<scalarField>
    (
        new scalarField(tkappa().boundaryField()[patchi])
    );
}


template<class BasicThermophysicalTransportModel>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::
legacyMixtureAverageFourier
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    Base(typeName, momentumTransport, thermo),
    TopoChangeableMeshObject(*this),
    DFuncs_(this->thermo().species().size()),
    legacyThermalDiffusionMode_("off"),
    Dm_(),
    sumRawFlux_()
{
    read();
}


template<class BasicThermophysicalTransportModel>
bool legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::read()
{
    if (!Base::read())
    {
        return false;
    }

    const speciesTable& species = this->thermo().species();
    const dictionary& coeffDict = this->typeDict();
    const dictionary& Ddict = coeffDict.subDict("D");

    legacyThermalDiffusionMode_ = coeffDict.lookupOrDefault<word>
    (
        "legacyThermalDiffusionMode",
        "off"
    );

    if (legacyThermalDiffusionMode_ != "off")
    {
        FatalIOErrorInFunction(coeffDict)
            << "D2 candidate-1 supports legacyThermalDiffusionMode off only. "
            << "The published OF8 H/H2 Soret implementation contains an "
            << "ambiguous TDRatio_H/TDRatio_H2 assignment and is being "
            << "qualified separately rather than silently changed."
            << exit(FatalIOError);
    }

    DFuncs_.setSize(species.size());

    forAll(species, i)
    {
        DFuncs_[i].setSize(species.size());

        for (label j = i + 1; j < species.size(); ++j)
        {
            const word nameij(species[i] + '-' + species[j]);
            const word nameji(species[j] + '-' + species[i]);
            word Dname;

            if (Ddict.found(nameij))
            {
                Dname = nameij;
            }
            else if (Ddict.found(nameji))
            {
                Dname = nameji;
            }
            else
            {
                FatalIOErrorInFunction(Ddict)
                    << "Missing binary diffusion coefficient for "
                    << nameij << " (or reversed name " << nameji << ')'
                    << exit(FatalIOError);
            }

            DFuncs_[i].set
            (
                j,
                Function2<scalar>::New
                (
                    Dname,
                    dimPressure,
                    dimTemperature,
                    dimKinematicViscosity,
                    Ddict
                ).ptr()
            );
        }
    }

    clearCache();
    return true;
}


template<class BasicThermophysicalTransportModel>
tmp<volScalarField>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::D
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "legacyD_" + Yi.name(),
        this->momentumTransport().rho()*Dm()[this->thermo().specieIndex(Yi)]
    );
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::D
(
    const volScalarField& Yi,
    const label patchi
) const
{
    const tmp<volScalarField> tDi = D(Yi);
    return tmp<scalarField>
    (
        new scalarField(tDi().boundaryField()[patchi])
    );
}


template<class BasicThermophysicalTransportModel>
tmp<volScalarField>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    return D(Yi);
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{
    return D(Yi, patchi);
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
    const label i = this->thermo().specieIndex(Yi);
    const tmp<surfaceScalarField> tRaw = rawFlux(i);
    const surfaceScalarField& sumRaw = sumRawFlux();

    return surfaceScalarField::New
    (
        IOobject::groupName("j(" + Yi.name() + ')', this->thermo().phaseName()),
       -tRaw()
      + fvc::interpolate(Yi)*sumRaw
    );
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::j
(
    const volScalarField& Yi,
    const label patchi
) const
{
    const tmp<surfaceScalarField> tJi = j(Yi);
    return tmp<scalarField>
    (
        new scalarField(tJi().boundaryField()[patchi])
    );
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    const volScalarField Wm(this->thermo().W());
    const tmp<volScalarField> tDi = D(Yi);
    const volScalarField& Di = tDi();

    // Keep the dominant Fickian term implicit, as in the OF8 equation.
    tmp<fvScalarMatrix> tDivj
    (
       -fvm::laplacian(this->alpha()*Di, Yi)
    );

    // Add the explicit molecular-weight and zero-net-mass-flux corrections.
    const tmp<volVectorField> tMw
    (
        Di*Yi/Wm*fvc::grad(Wm)
    );

    const surfaceScalarField mwFlux
    (
        "legacyMwFlux_" + Yi.name(),
        (
            fvc::interpolate(this->alpha()*tMw())
          & Yi.mesh().Sf()
        )/Yi.mesh().magSf()
    );

    const surfaceScalarField& sumRaw = sumRawFlux();
    const surfaceScalarField jCorrection
    (
        "legacyJCorrection_" + Yi.name(),
       -mwFlux + fvc::interpolate(Yi)*sumRaw
    );

    tDivj.ref() += fvc::div(jCorrection*Yi.mesh().magSf());
    return tDivj;
}


template<class BasicThermophysicalTransportModel>
tmp<surfaceScalarField>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::q() const
{
    tmp<surfaceScalarField> tq
    (
        surfaceScalarField::New
        (
            IOobject::groupName("q", this->thermo().phaseName()),
           -fvc::interpolate(this->alpha()*this->legacyKappa())
           *fvc::snGrad(this->thermo().T())
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();

    surfaceScalarField sumJ
    (
        surfaceScalarField::New
        (
            "legacySumJ",
            T.mesh(),
            dimensionedScalar(dimMass/dimArea/dimTime, 0)
        )
    );
    surfaceScalarField sumJh
    (
        surfaceScalarField::New
        (
            "legacySumJh",
            T.mesh(),
            dimensionedScalar(dimEnergy/dimArea/dimTime, 0)
        )
    );

    forAll(Y, i)
    {
        if (i != this->thermo().defaultSpecie())
        {
            const volScalarField hi(this->thermo().hsi(i, p, T));
            const surfaceScalarField ji(this->j(Y[i]));
            sumJ += ji;
            sumJh += ji*fvc::interpolate(hi);
        }
    }

    const label di = this->thermo().defaultSpecie();
    const volScalarField hDefault(this->thermo().hsi(di, p, T));
    sumJh -= sumJ*fvc::interpolate(hDefault);
    tq.ref() += sumJh;

    return tq;
}


template<class BasicThermophysicalTransportModel>
tmp<scalarField>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::q
(
    const label patchi
) const
{
    // Match the native OF14 Fickian patch heat-flux expression exactly in
    // operator grouping.  alpha() is oneField for this instantiation, so the
    // unary minus must apply to the complete scalar-field product rather than
    // directly to oneField.
    tmp<scalarField> tq
    (
       -(
            this->alpha().boundaryField()[patchi]
           *this->legacyKappa(patchi)
           *this->thermo().T().boundaryField()[patchi].snGrad()
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const scalarField& p = this->thermo().p().boundaryField()[patchi];
    const scalarField& T = this->thermo().T().boundaryField()[patchi];
    scalarField sumJ(tq().size(), scalar(0));
    scalarField sumJh(tq().size(), scalar(0));

    forAll(Y, i)
    {
        if (i != this->thermo().defaultSpecie())
        {
            const scalarField hi(this->thermo().hsi(i, p, T));
            const tmp<scalarField> tji = j(Y[i], patchi);
            sumJ += tji();
            sumJh += tji()*hi;
        }
    }

    const label di = this->thermo().defaultSpecie();
    sumJh -= sumJ*this->thermo().hsi(di, p, T);
    tq.ref() += sumJh;
    return tq;
}


template<class BasicThermophysicalTransportModel>
tmp<fvScalarMatrix>
legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
    const tmp<volScalarField> tkappa(this->legacyKappa());
    const volScalarField& kappa = tkappa();

    tmp<fvScalarMatrix> tDivq
    (
        fvm::Su
        (
           -fvc::laplacian
            (
                this->alpha()*kappa,
                this->thermo().T()
            ),
            he
        )
    );

    tDivq.ref() -=
        fvm::laplacianCorrection
        (
            this->alpha()*kappa/this->thermo().Cpv(),
            he
        );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();
    surfaceScalarField sumJ
    (
        surfaceScalarField::New
        (
            "legacySumJ",
            he.mesh(),
            dimensionedScalar(dimMass/dimArea/dimTime, 0)
        )
    );
    surfaceScalarField sumJh
    (
        surfaceScalarField::New
        (
            "legacySumJh",
            he.mesh(),
            dimensionedScalar(dimEnergy/dimArea/dimTime, 0)
        )
    );

    forAll(Y, i)
    {
        if (i != this->thermo().defaultSpecie())
        {
            const volScalarField hi(this->thermo().hsi(i, p, T));
            const surfaceScalarField ji(this->j(Y[i]));
            sumJ += ji;
            sumJh += ji*fvc::interpolate(hi);
        }
    }

    const label di = this->thermo().defaultSpecie();
    const volScalarField hDefault(this->thermo().hsi(di, p, T));
    sumJh -= sumJ*fvc::interpolate(hDefault);

    tDivq.ref() += fvc::div(sumJh*he.mesh().magSf());
    return tDivq;
}


template<class BasicThermophysicalTransportModel>
void legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::predict()
{
    Base::predict();
    clearCache();
    updateDm();
    updateSumRawFlux();
}


template<class BasicThermophysicalTransportModel>
bool legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::movePoints()
{
    clearCache();
    return true;
}


template<class BasicThermophysicalTransportModel>
void legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::topoChange
(
    const polyTopoChangeMap&
)
{
    clearCache();
}


template<class BasicThermophysicalTransportModel>
void legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::mapMesh
(
    const polyMeshMap&
)
{
    clearCache();
}


template<class BasicThermophysicalTransportModel>
void legacyMixtureAverageFourier<BasicThermophysicalTransportModel>::distribute
(
    const polyDistributionMap&
)
{
    clearCache();
}

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam
