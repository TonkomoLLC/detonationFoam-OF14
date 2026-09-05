#include "planarRefiner.H"
#include "planarMultiDirRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "polyBoundaryMesh.H"
#include "DynamicList.H"
#include "PstreamReduceOps.H"

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(planarRefiner, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, planarRefiner, fvMesh);
}
}

Foam::word Foam::fvMeshTopoChangers::planarRefiner::detectedGeometry() const
{
    label nEmpty=0;
    label nWedge=0;

    const polyBoundaryMesh& patches=mesh().poly().boundary();
    forAll(patches, patchi)
    {
        const word patchType=patches[patchi].type();
        if (patchType == "empty") ++nEmpty;
        if (patchType == "wedge") ++nWedge;
    }

    if (nWedge >= 2 && nEmpty == 0) return "wedge";
    if (nEmpty >= 1 && nWedge == 0) return "slab";

    FatalErrorInFunction
        << "Cannot identify a supported two-dimensional geometry from boundary patches. "
        << "Found " << nEmpty << " empty and " << nWedge << " wedge patches. "
        << "Supported geometries are slab/empty and wedge/axisymmetric."
        << exit(FatalError);

    return word::null;
}

void Foam::fvMeshTopoChangers::planarRefiner::validateGeometry()
{
    if (mesh().nGeometricD() != 2)
    {
        FatalIOErrorInFunction(dict_)
            << "planarRefiner requires a 2-D mesh; nGeometricD="
            << mesh().nGeometricD() << exit(FatalIOError);
    }

    if (geometry_ != "auto" && geometry_ != "slab" && geometry_ != "wedge")
    {
        FatalIOErrorInFunction(dict_)
            << "geometry must be auto, slab, or wedge; got " << geometry_
            << exit(FatalIOError);
    }

    const word detected=detectedGeometry();
    if (geometry_ == "auto")
    {
        geometry_=detected;
    }
    else if (geometry_ != detected)
    {
        FatalIOErrorInFunction(dict_)
            << "Requested geometry " << geometry_
            << " does not match boundary-patch geometry " << detected
            << exit(FatalIOError);
    }
}

Foam::dictionary Foam::fvMeshTopoChangers::planarRefiner::directionalRefineDict() const
{
    if (mesh().nGeometricD() != 2)
    {
        FatalErrorInFunction
            << "planarRefiner requires a two-dimensional fvMesh; nGeometricD="
            << mesh().nGeometricD() << exit(FatalError);
    }

    // OpenFOAM's geometricD() identifies the inactive direction for both
    // empty/slab and wedge/axisymmetric meshes.  Refine the two active
    // directions only; the wedge/empty constraint is corrected by fvMesh's
    // normal topology-change machinery after each directional split.
    const Vector<label> dirs(mesh().geometricD());
    wordList directions(2);
    if (dirs.x() == -1)
    {
        directions[0]="e2";
        directions[1]="e3";
    }
    else if (dirs.y() == -1)
    {
        directions[0]="e1";
        directions[1]="e3";
    }
    else
    {
        directions[0]="e1";
        directions[1]="e2";
    }

    dictionary d;
    d.add("directions", directions);
    d.add("useHexTopology", false);
    d.add("coordinateSystem", word("global"));
    dictionary global;
    global.add("e1", vector(1,0,0));
    global.add("e2", vector(0,1,0));
    d.add("global", global);
    d.add("geometricCut", false);
    d.add("writeMesh", false);
    return d;
}

Foam::labelList Foam::fvMeshTopoChangers::planarRefiner::selectCells() const
{
    if (!mesh().foundObject<volScalarField>(fieldName_))
    {
        FatalErrorInFunction
            << "Refinement field " << fieldName_
            << " is not registered as a volScalarField"
            << exit(FatalError);
    }

    const volScalarField& fld=mesh().lookupObject<volScalarField>(fieldName_);
    DynamicList<label> selected;
    forAll(fld,celli)
    {
        if (fld[celli] >= lowerRefineLevel_ && fld[celli] <= upperRefineLevel_)
        {
            selected.append(celli);
        }
    }
    return labelList(selected);
}

Foam::fvMeshTopoChangers::planarRefiner::planarRefiner
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshTopoChanger(mesh),
    dict_(dict),
    refineInterval_(dict.lookup<label>("refineInterval")),
    fieldName_(dict.lookup<word>("field")),
    lowerRefineLevel_(dict.lookup<scalar>("lowerRefineLevel")),
    upperRefineLevel_(dict.lookup<scalar>("upperRefineLevel")),
    maxCells_(dict.lookup<label>("maxCells")),
    maxRefinementIterations_(dict.lookupOrDefault<label>("maxRefinementIterations",1)),
    nRefinementIterations_(0),
    geometry_(dict.lookupOrDefault<word>("geometry","auto")),
    timeIndex_(-1),
    state_
    (
        IOobject
        (
            "planarRefinerState",
            mesh.time().name(),
            "polyMesh",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    )
{
    if (refineInterval_ < 1)
    {
        FatalIOErrorInFunction(dict_)
            << "refineInterval must be >= 1" << exit(FatalIOError);
    }
    if (maxCells_ <= 0 || maxRefinementIterations_ < 1)
    {
        FatalIOErrorInFunction(dict_)
            << "maxCells must be > 0 and maxRefinementIterations >= 1"
            << exit(FatalIOError);
    }

    validateGeometry();

    nRefinementIterations_=
        state_.lookupOrDefault<label>("nRefinementIterations",0);

    label minIter=nRefinementIterations_;
    label maxIter=nRefinementIterations_;
    reduce(minIter,minOp<label>());
    reduce(maxIter,maxOp<label>());
    if (minIter != maxIter)
    {
        FatalErrorInFunction
            << "Inconsistent planarRefiner restart state across MPI ranks: min="
            << minIter << " max=" << maxIter << exit(FatalError);
    }
    if (nRefinementIterations_ < 0
     || nRefinementIterations_ > maxRefinementIterations_)
    {
        FatalErrorInFunction
            << "Invalid persisted nRefinementIterations="
            << nRefinementIterations_ << "; configured maximum is "
            << maxRefinementIterations_ << exit(FatalError);
    }

    Info<< "planarRefiner: reusable OF14 2-D slab/wedge dynamic refiner" << nl
        << "    geometry=" << geometry_ << nl
        << "    field=" << fieldName_ << " band="
        << lowerRefineLevel_ << ".." << upperRefineLevel_ << nl
        << "    maxCells=" << maxCells_
        << " maxRefinementIterations=" << maxRefinementIterations_ << nl
        << "    restored refinement iteration " << nRefinementIterations_
        << "/" << maxRefinementIterations_ << endl;
}

Foam::fvMeshTopoChangers::planarRefiner::~planarRefiner() {}

bool Foam::fvMeshTopoChangers::planarRefiner::update()
{
    const label ti=mesh().time().timeIndex();
    if (timeIndex_ == ti) return false;
    timeIndex_=ti;

    if
    (
        ti == 0
     || (ti % refineInterval_) != 0
     || nRefinementIterations_ >= maxRefinementIterations_
    )
    {
        return false;
    }

    const label oldCells=mesh().globalData().nTotalCells();
    if (oldCells >= maxCells_) return false;

    labelList cellsToRefine(selectCells());
    const label nSelected=returnReduce(cellsToRefine.size(),sumOp<label>());

    if (!nSelected)
    {
        Info<< "planarRefiner: no cells selected at time "
            << mesh().time().name() << endl;
        return false;
    }

    // A planar one-level split creates four children per selected cell,
    // i.e. a net +3 cells for each selected parent.  The directional cutter
    // does not add hidden refinement levels, so this is an exact growth bound
    // for the current one-pass implementation.
    if (oldCells + 3*nSelected > maxCells_)
    {
        Info<< "planarRefiner: planar growth would exceed maxCells; "
            << "refinement skipped" << endl;
        return false;
    }

    dictionary refineDict(directionalRefineDict());
    Info<< "planarRefiner: refining " << nSelected
        << " localized cells at time " << mesh().time().name() << endl;

    planarMultiDirRefinement(mesh(),cellsToRefine,refineDict,refineDict);

    const label newCells=mesh().globalData().nTotalCells();
    if (newCells != oldCells)
    {
        ++nRefinementIterations_;
    }

    Info<< "Planar refined from " << oldCells << " to " << newCells
        << " cells." << nl
        << "planarRefiner: refinement iteration " << nRefinementIterations_
        << "/" << maxRefinementIterations_ << endl;

    return newCells != oldCells;
}

void Foam::fvMeshTopoChangers::planarRefiner::topoChange
(
    const polyTopoChangeMap&
)
{
    // No independent per-cell refinement metadata is stored here.  Registered
    // fvMesh fields are mapped by fvMesh::topoChange() after each directional
    // split in planarRefinementIterator.
}

void Foam::fvMeshTopoChangers::planarRefiner::mapMesh(const polyMeshMap&)
{
    // The only persistent state is the global refinement-iteration counter.
}

void Foam::fvMeshTopoChangers::planarRefiner::distribute
(
    const polyDistributionMap&
)
{
    // Native fvMesh distribution owns rank-local mesh/field redistribution.
    // planarRefiner has no cell-indexed private state to distribute.
    Info<< "planarRefiner: received native mesh-distribution callback" << endl;
}

bool Foam::fvMeshTopoChangers::planarRefiner::write(const bool doWrite) const
{
    if (!doWrite) return true;

    state_.set("formatVersion",1);
    state_.set("nRefinementIterations",nRefinementIterations_);
    state_.set("globalCells",mesh().globalData().nTotalCells());
    state_.set("geometry",geometry_);

    // The iteration counter is global (the same on every rank), and
    // IOdictionary is a global OpenFOAM IO object.  Write one case-level
    // state dictionary at the current time while OpenFOAM owns the rank-local
    // distributed mesh and registered field state.
    state_.instance()=mesh().time().name();
    const bool ok=state_.regIOobject::write(doWrite);

    if (ok)
    {
        Info<< "planarRefiner: wrote restart state at time "
            << mesh().time().name() << " iteration="
            << nRefinementIterations_ << endl;
    }
    return ok;
}
