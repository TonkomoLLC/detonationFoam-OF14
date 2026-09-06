/*---------------------------------------------------------------------------*\
  Qualification-only driver for planarRefiner automatic unrefinement.
  It registers deterministic scalar fields, changes only the refinement
  indicator, calls fvMesh::update(), and reports global mapping invariants.
\*---------------------------------------------------------------------------*/
#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "PstreamReduceOps.H"

using namespace Foam;

namespace
{
struct Metrics
{
    label cells;
    scalar volume;
    scalar uniformIntegral;
    scalar tracerIntegral;
    scalar uniformMin;
    scalar uniformMax;
};

Metrics metrics
(
    const fvMesh& mesh,
    const volScalarField& uniformTracer,
    const volScalarField& tracer
)
{
    scalar v = 0;
    scalar ui = 0;
    scalar ti = 0;
    scalar uMin = GREAT;
    scalar uMax = -GREAT;

    const auto& V = mesh.V();
    forAll(V, celli)
    {
        v += V[celli];
        ui += V[celli]*uniformTracer[celli];
        ti += V[celli]*tracer[celli];
        uMin = min(uMin, uniformTracer[celli]);
        uMax = max(uMax, uniformTracer[celli]);
    }

    reduce(v, sumOp<scalar>());
    reduce(ui, sumOp<scalar>());
    reduce(ti, sumOp<scalar>());
    reduce(uMin, minOp<scalar>());
    reduce(uMax, maxOp<scalar>());

    Metrics m;
    m.cells = mesh.globalData().nTotalCells();
    m.volume = v;
    m.uniformIntegral = ui;
    m.tracerIntegral = ti;
    m.uniformMin = uMin;
    m.uniformMax = uMax;
    return m;
}

void printMetrics(const word& tag, const label step, const Metrics& m)
{
    Info<< "QUAL_METRIC tag=" << tag
        << " step=" << step
        << " cells=" << m.cells
        << " volume=" << m.volume
        << " uniformIntegral=" << m.uniformIntegral
        << " tracerIntegral=" << m.tracerIntegral
        << " uniformMin=" << m.uniformMin
        << " uniformMax=" << m.uniformMax
        << endl;
}
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary qDict
    (
        IOobject
        (
            "qualificationProperties",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word mode(qDict.lookupOrDefault<word>("mode", "uniform"));
    const scalar refineValue=qDict.lookupOrDefault<scalar>("refineValue",1.0);
    const scalar unrefineValue=qDict.lookupOrDefault<scalar>("unrefineValue",0.0);
    const scalar xMin=qDict.lookupOrDefault<scalar>("xMin",-GREAT);
    const scalar xMax=qDict.lookupOrDefault<scalar>("xMax", GREAT);
    const scalar uniformValue=qDict.lookupOrDefault<scalar>("uniformTracerValue",2.5);

    volScalarField indicator
    (
        IOobject
        (
            "indicator",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField uniformTracer
    (
        IOobject
        (
            "uniformTracer",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField tracer
    (
        IOobject
        (
            "tracer",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Only initialise the mapping probes for a base-time run. A restart reads
    // the already mapped fields and must not overwrite them.
    if (runTime.value() == 0)
    {
        const volVectorField& C = mesh.C();
        forAll(tracer, celli)
        {
            uniformTracer[celli] = uniformValue;
            tracer[celli] = 1.0 + C[celli].x() + 0.5*C[celli].y();
            indicator[celli] = unrefineValue;
        }
        uniformTracer.correctBoundaryConditions();
        tracer.correctBoundaryConditions();
        indicator.correctBoundaryConditions();
    }

    const Metrics base=metrics(mesh,uniformTracer,tracer);
    printMetrics("base",runTime.timeIndex(),base);

    while (runTime.run())
    {
        ++runTime;
        const label step=runTime.timeIndex();
        const bool refinePhase=(step % 2) == 1;
        const volVectorField& C = mesh.C();

        forAll(indicator, celli)
        {
            bool selected=true;
            if (mode == "localized")
            {
                selected = C[celli].x() >= xMin && C[celli].x() <= xMax;
            }
            else if (mode != "uniform")
            {
                FatalIOErrorInFunction(qDict)
                    << "mode must be uniform or localized; got " << mode
                    << exit(FatalIOError);
            }

            indicator[celli] =
                (refinePhase && selected) ? refineValue : unrefineValue;
        }
        indicator.correctBoundaryConditions();

        const Metrics before=metrics(mesh,uniformTracer,tracer);
        printMetrics(refinePhase ? "beforeRefine" : "beforeUnrefine",step,before);

        const bool changed=mesh.update();

        const Metrics after=metrics(mesh,uniformTracer,tracer);
        printMetrics(refinePhase ? "afterRefine" : "afterUnrefine",step,after);
        Info<< "QUAL_CHANGE step=" << step
            << " changed=" << changed
            << " cellsBefore=" << before.cells
            << " cellsAfter=" << after.cells
            << endl;

        runTime.write();
    }

    Info<< "QUAL_END" << endl;
    return 0;
}
