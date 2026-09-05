/*---------------------------------------------------------------------------*\
  legacyMixtureAverageFourier
  Minimal OpenFOAM Foundation 14 transport model for detonationFoam legacy
  mixture-averaged species diffusion compatibility.
\*---------------------------------------------------------------------------*/
#ifndef legacyMixtureAverageFourier_H
#define legacyMixtureAverageFourier_H

#include "unityLewisFourier.H"
#include "Function2.H"

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

template<class BasicThermophysicalTransportModel>
class legacyMixtureAverageFourier
:
    public unityLewisFourier<BasicThermophysicalTransportModel>,
    public TopoChangeableMeshObject<fvMesh>
{
    typedef unityLewisFourier<BasicThermophysicalTransportModel> Base;

    // Binary D_ij(p,T) functions. Only the upper triangle i<j is populated.
    List<PtrList<Function2<scalar>>> DFuncs_;

    // D2 candidate-1 deliberately supports only the historically exercised
    // thermalDiffusion=off path. The published OF8 H2 Soret assignment is
    // ambiguous and is not silently corrected here.
    word legacyThermalDiffusionMode_;

    // Per-predictor caches. They are cleared on every mesh topology/map/
    // distribution change, mirroring the native OF14 Fickian cache policy.
    mutable PtrList<volScalarField> Dm_;
    mutable tmp<surfaceScalarField> sumRawFlux_;

    void clearCache() const;
    void updateDm() const;
    const PtrList<volScalarField>& Dm() const;

    const Function2<scalar>& Dij(const label i, const label j) const;

    // Positive legacy raw diffusion vector before the mixture correction:
    // rho*Di*grad(Yi) + Yi*rho*Di*grad(W)/W.
    tmp<surfaceScalarField> rawFlux(const label i) const;

    void updateSumRawFlux() const;
    const surfaceScalarField& sumRawFlux() const;

    // Exact legacy OF8 thermal-conductivity mixture rule:
    // 0.5*(sum_i Xi*kappa_i + 1/sum_i(Xi/kappa_i)).
    // Species kappa_i are provided by the selected OF14 thermo package.
    tmp<volScalarField> legacyKappa() const;
    tmp<scalarField> legacyKappa(const label patchi) const;

public:
    typedef typename BasicThermophysicalTransportModel::alphaField alphaField;
    typedef typename BasicThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;
    typedef typename BasicThermophysicalTransportModel::thermoModel thermoModel;

    TypeName("legacyMixtureAverageFourier");

    legacyMixtureAverageFourier
    (
        const momentumTransportModel& momentumTransport,
        const thermoModel& thermo
    );

    virtual ~legacyMixtureAverageFourier() {}

    virtual bool read();

    virtual tmp<volScalarField> D(const volScalarField& Yi) const;
    virtual tmp<scalarField> D(const volScalarField& Yi, const label patchi) const;

    virtual tmp<volScalarField> DEff(const volScalarField& Yi) const;
    virtual tmp<scalarField> DEff(const volScalarField& Yi, const label patchi) const;

    virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const;
    virtual tmp<scalarField> j(const volScalarField& Yi, const label patchi) const;
    virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

    virtual tmp<surfaceScalarField> q() const;
    virtual tmp<scalarField> q(const label patchi) const;
    virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

    virtual void predict();

    virtual bool movePoints();
    virtual void topoChange(const polyTopoChangeMap& map);
    virtual void mapMesh(const polyMeshMap& map);
    virtual void distribute(const polyDistributionMap& map);
};

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

#ifdef NoRepository
    #include "legacyMixtureAverageFourier.C"
#endif

#endif
