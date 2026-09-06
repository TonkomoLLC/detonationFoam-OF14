#include "planarRefiner.H"
#include "planarMultiDirRefinement.H"
#include "undoableMeshCutter.H"
#include "polyTopoChange.H"
#include "polyTopoChangeMap.H"
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


void Foam::fvMeshTopoChangers::planarRefiner::validateAutomaticUnrefinement()
{
    if (!automaticUnrefinement_)
    {
        return;
    }

    if (unrefineInterval_ < 1)
    {
        FatalIOErrorInFunction(dict_)
            << "unrefineInterval must be >= 1" << exit(FatalIOError);
    }

    if (maxUnrefinementPassesPerUpdate_ < 1)
    {
        FatalIOErrorInFunction(dict_)
            << "maxUnrefinementPassesPerUpdate must be >= 1"
            << exit(FatalIOError);
    }

    if (maxRefinementLevel_ < 1)
    {
        FatalIOErrorInFunction(dict_)
            << "maxRefinementLevel must be >= 1" << exit(FatalIOError);
    }

    if (lowerUnrefineLevel_ > lowerRefineLevel_)
    {
        FatalIOErrorInFunction(dict_)
            << "lowerUnrefineLevel=" << lowerUnrefineLevel_
            << " must be <= lowerRefineLevel=" << lowerRefineLevel_
            << ". This ordering provides the low-side hysteresis band."
            << exit(FatalIOError);
    }

    if (upperUnrefineLevel_ < upperRefineLevel_)
    {
        FatalIOErrorInFunction(dict_)
            << "upperUnrefineLevel=" << upperUnrefineLevel_
            << " must be >= upperRefineLevel=" << upperRefineLevel_
            << ". This ordering provides the high-side hysteresis band."
            << exit(FatalIOError);
    }

    if (lowerUnrefineLevel_ >= upperUnrefineLevel_)
    {
        FatalIOErrorInFunction(dict_)
            << "lowerUnrefineLevel must be smaller than upperUnrefineLevel"
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
    // empty/slab and wedge/axisymmetric meshes. Refine the two active
    // directions only.
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
        const bool inBand =
            fld[celli] >= lowerRefineLevel_
         && fld[celli] <= upperRefineLevel_;

        // One planar refinement level is two binary directional cuts.
        // Only select a cell when a complete additional planar level fits
        // below the requested depth cap.  This prevents a partial previous
        // cut from permitting a later two-cut operation to overshoot it.
        const bool belowLevelCap =
            !automaticUnrefinement_
         || refinementLevel_[celli] + 2
                <= 2*maxRefinementLevel_ + SMALL;

        if (inBand && belowLevelCap)
        {
            selected.append(celli);
        }
    }
    return labelList(selected);
}


Foam::labelList
Foam::fvMeshTopoChangers::planarRefiner::selectSplitFacesForUnrefinement() const
{
    if (!automaticUnrefinement_ || !reversibleCutter_.valid())
    {
        return labelList();
    }

    if (!mesh().foundObject<volScalarField>(fieldName_))
    {
        FatalErrorInFunction
            << "Unrefinement field " << fieldName_
            << " is not registered as a volScalarField"
            << exit(FatalError);
    }

    const volScalarField& fld=mesh().lookupObject<volScalarField>(fieldName_);
    const labelList splitFaces(reversibleCutter_->getSplitFaces());
    const labelUList& owner=mesh().faceOwner();
    const labelUList& neighbour=mesh().faceNeighbour();

    DynamicList<label> selected(splitFaces.size());

    forAll(splitFaces, i)
    {
        const label facei=splitFaces[i];

        // A reversible split face produced inside a parent cell must be an
        // internal face. Refuse an unexpected topology instead of guessing.
        if (facei < 0 || facei >= mesh().nInternalFaces())
        {
            FatalErrorInFunction
                << "undoableMeshCutter returned non-internal split face "
                << facei << ". Automatic planar unrefinement is refusing "
                << "this topology." << exit(FatalError);
        }

        const label own=owner[facei];
        const label nei=neighbour[facei];

        // Only merge siblings at the same non-zero binary cut depth.
        if
        (
            refinementLevel_[own] < 0.5
         || refinementLevel_[nei] < 0.5
         || mag(refinementLevel_[own] - refinementLevel_[nei]) > SMALL
        )
        {
            continue;
        }

        // Conservative two-sided hysteresis: both siblings must have left
        // the active refinement band on the same side. This prevents merging
        // a face that still straddles a shock/reaction-front indicator.
        const bool lowSide =
            fld[own] < lowerUnrefineLevel_
         && fld[nei] < lowerUnrefineLevel_;
        const bool highSide =
            fld[own] > upperUnrefineLevel_
         && fld[nei] > upperUnrefineLevel_;

        if (lowSide || highSide)
        {
            selected.append(facei);
        }
    }

    return labelList(selected);
}


bool Foam::fvMeshTopoChangers::planarRefiner::unrefineOnePass()
{
    labelList splitFaces(selectSplitFacesForUnrefinement());
    const label nSelected=returnReduce(splitFaces.size(),sumOp<label>());

    if (!nSelected)
    {
        return false;
    }

    const scalarField oldLevels(refinementLevel_.primitiveField());
    const labelUList& owner=mesh().faceOwner();
    const labelUList& neighbour=mesh().faceNeighbour();

    labelList affectedCells(2*splitFaces.size());
    forAll(splitFaces, i)
    {
        affectedCells[2*i]=owner[splitFaces[i]];
        affectedCells[2*i+1]=neighbour[splitFaces[i]];
    }

    const label oldCells=mesh().globalData().nTotalCells();

    mesh().preChange();
    polyTopoChange meshMod(mesh());

    const labelList removedFaces
    (
        reversibleCutter_->removeSplitFaces(splitFaces,meshMod)
    );

    // OpenFOAM documents that removeFaces can occasionally request extra
    // faces. That case is deliberately rejected here because it breaks the
    // one-split-face/one-binary-level bookkeeping used by planarRefiner.
    if (removedFaces.size() != splitFaces.size())
    {
        FatalErrorInFunction
            << "Automatic planar unrefinement requested " << splitFaces.size()
            << " split faces locally, but undoableMeshCutter/removeFaces "
            << "selected " << removedFaces.size() << ". Version 1.1.0 refuses "
            << "additional face removal because refinement-level bookkeeping "
            << "would no longer be exact." << exit(FatalError);
    }

    autoPtr<polyTopoChangeMap> map=meshMod.changeMesh(mesh());

    mesh().topoChange(map());
    reversibleCutter_->topoChange(map());

    const labelList& reverseCellMap=map->reverseCellMap();
    forAll(affectedCells, i)
    {
        const label oldCelli=affectedCells[i];
        if (oldCelli >= 0 && oldCelli < reverseCellMap.size())
        {
            const label newCelli=reverseCellMap[oldCelli];
            if (newCelli >= 0 && newCelli < refinementLevel_.size())
            {
                refinementLevel_[newCelli]=max(scalar(0),oldLevels[oldCelli]-1);
            }
        }
    }
    refinementLevel_.correctBoundaryConditions();

    const label newCells=mesh().globalData().nTotalCells();
    if (newCells < oldCells)
    {
        ++nUnrefinementPasses_;
        Info<< "planarRefiner: unrefined " << nSelected
            << " reversible split faces; global cells " << oldCells
            << " -> " << newCells << endl;
        return true;
    }

    FatalErrorInFunction
        << "undoableMeshCutter accepted " << nSelected
        << " split faces but global cell count did not decrease ("
        << oldCells << " -> " << newCells << ')'
        << exit(FatalError);

    return false;
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
    automaticUnrefinement_
    (
        dict.lookupOrDefault<Switch>("automaticUnrefinement",false)
    ),
    unrefineInterval_
    (
        dict.lookupOrDefault<label>("unrefineInterval",refineInterval_)
    ),
    lowerUnrefineLevel_
    (
        dict.lookupOrDefault<scalar>("lowerUnrefineLevel",lowerRefineLevel_)
    ),
    upperUnrefineLevel_
    (
        dict.lookupOrDefault<scalar>("upperUnrefineLevel",upperRefineLevel_)
    ),
    maxUnrefinementPassesPerUpdate_
    (
        dict.lookupOrDefault<label>("maxUnrefinementPassesPerUpdate",2)
    ),
    maxRefinementLevel_
    (
        dict.lookupOrDefault<label>
        (
            "maxRefinementLevel",
            maxRefinementIterations_
        )
    ),
    nUnrefinementPasses_(0),
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
    ),
    reversibleCutter_(nullptr),
    refinementLevel_
    (
        IOobject
        (
            "planarRefinementLevel",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimless,0)
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
    if (lowerRefineLevel_ >= upperRefineLevel_)
    {
        FatalIOErrorInFunction(dict_)
            << "lowerRefineLevel must be smaller than upperRefineLevel"
            << exit(FatalIOError);
    }

    validateGeometry();
    validateAutomaticUnrefinement();

    nRefinementIterations_=
        state_.lookupOrDefault<label>("nRefinementIterations",0);
    nUnrefinementPasses_=
        state_.lookupOrDefault<label>("nUnrefinementPasses",0);

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

    if (nRefinementIterations_ < 0)
    {
        FatalErrorInFunction
            << "Invalid persisted nRefinementIterations="
            << nRefinementIterations_ << exit(FatalError);
    }

    if
    (
        !automaticUnrefinement_
     && nRefinementIterations_ > maxRefinementIterations_
    )
    {
        FatalErrorInFunction
            << "Invalid persisted nRefinementIterations="
            << nRefinementIterations_ << "; configured maximum is "
            << maxRefinementIterations_ << exit(FatalError);
    }

    if (automaticUnrefinement_)
    {
        const label formatVersion=state_.lookupOrDefault<label>("formatVersion",1);
        const bool activeHistory=state_.lookupOrDefault<bool>
        (
            "activeReversibleSplits",
            formatVersion < 2 && nRefinementIterations_ > 0
        );

        if (activeHistory)
        {
            FatalErrorInFunction
                << "automaticUnrefinement is enabled, but this restart state "
                << "contains active reversible planar refinement. Version 1.1.0 "
                << "does not serialize undoableMeshCutter history, so an exact "
                << "restart cannot reconstruct safe unrefinement ancestry. "
                << "Restart from a time written after all reversible cuts have "
                << "been coarsened, or set automaticUnrefinement false."
                << exit(FatalError);
        }

        reversibleCutter_.reset(new undoableMeshCutter(mesh,true));
    }

    Info<< "planarRefiner: reusable OF14 2-D slab/wedge dynamic refiner" << nl
        << "    geometry=" << geometry_ << nl
        << "    field=" << fieldName_ << " refine band="
        << lowerRefineLevel_ << ".." << upperRefineLevel_ << nl
        << "    maxCells=" << maxCells_;

    if (automaticUnrefinement_)
    {
        Info<< nl
            << "    automaticUnrefinement=true, unrefine outside "
            << lowerUnrefineLevel_ << ".." << upperUnrefineLevel_ << nl
            << "    maxRefinementLevel=" << maxRefinementLevel_
            << " maxUnrefinementPassesPerUpdate="
            << maxUnrefinementPassesPerUpdate_ << nl
            << "    NOTE: active reversible-history restart and mesh "
            << "redistribution/load balancing are disabled in v1.1.0";
    }
    else
    {
        Info<< " maxRefinementIterations=" << maxRefinementIterations_ << nl
            << "    restored refinement iteration " << nRefinementIterations_
            << "/" << maxRefinementIterations_;
    }
    Info<< endl;
}


Foam::fvMeshTopoChangers::planarRefiner::~planarRefiner() {}


bool Foam::fvMeshTopoChangers::planarRefiner::update()
{
    const label ti=mesh().time().timeIndex();
    if (timeIndex_ == ti) return false;
    timeIndex_=ti;

    if (ti == 0)
    {
        return false;
    }

    bool changed=false;

    // Coarsen first, then permit refinement. For a moving front this frees
    // cells behind/ahead of the active indicator before maxCells is checked.
    if (automaticUnrefinement_ && (ti % unrefineInterval_) == 0)
    {
        for (label pass=0; pass<maxUnrefinementPassesPerUpdate_; ++pass)
        {
            if (!unrefineOnePass())
            {
                break;
            }
            changed=true;
        }
    }

    if ((ti % refineInterval_) != 0)
    {
        return changed;
    }

    // Preserve the original global-iteration behavior exactly when automatic
    // unrefinement is disabled. Reversible AMR instead uses per-cell depth.
    if
    (
        !automaticUnrefinement_
     && nRefinementIterations_ >= maxRefinementIterations_
    )
    {
        return changed;
    }

    const label oldCells=mesh().globalData().nTotalCells();
    if (oldCells >= maxCells_)
    {
        return changed;
    }

    labelList cellsToRefine(selectCells());
    const label nSelected=returnReduce(cellsToRefine.size(),sumOp<label>());

    if (!nSelected)
    {
        Info<< "planarRefiner: no cells selected for refinement at time "
            << mesh().time().name() << endl;
        return changed;
    }

    // A complete planar one-level split creates four children per selected
    // cell, i.e. net +3. This is the exact maximum for the two directional
    // cuts used here.
    if (oldCells + 3*nSelected > maxCells_)
    {
        Info<< "planarRefiner: planar growth would exceed maxCells; "
            << "refinement skipped" << endl;
        return changed;
    }

    dictionary refineDict(directionalRefineDict());
    Info<< "planarRefiner: refining " << nSelected
        << " localized cells at time " << mesh().time().name() << endl;

    const scalarField oldLevels(refinementLevel_.primitiveField());

    if (automaticUnrefinement_)
    {
        planarMultiDirRefinement refiner
        (
            mesh(),
            reversibleCutter_(),
            cellsToRefine,
            refineDict,
            refineDict
        );

        // Registered fields have already been mapped by fvMesh::topoChange()
        // after each directional cut. Set exact binary depth on all final
        // children originating from the selected old parents.
        const labelListList& added=refiner.addedCells();
        forAll(cellsToRefine, i)
        {
            const label oldCelli=cellsToRefine[i];
            if (oldCelli < 0 || oldCelli >= added.size()) continue;

            const labelList& children=added[oldCelli];
            if (children.empty()) continue;

            // Two children = one successful binary directional cut;
            // four children = both planar directions succeeded.
            scalar increment=0;
            if (children.size() >= 4) increment=2;
            else if (children.size() >= 2) increment=1;

            forAll(children, childI)
            {
                const label celli=children[childI];
                if (celli >= 0 && celli < refinementLevel_.size())
                {
                    refinementLevel_[celli]=oldLevels[oldCelli]+increment;
                }
            }
        }
        refinementLevel_.correctBoundaryConditions();
    }
    else
    {
        planarMultiDirRefinement(mesh(),cellsToRefine,refineDict,refineDict);
    }

    const label newCells=mesh().globalData().nTotalCells();
    if (newCells != oldCells)
    {
        ++nRefinementIterations_;
        changed=true;
    }

    Info<< "Planar refined from " << oldCells << " to " << newCells
        << " cells." << nl;
    if (automaticUnrefinement_)
    {
        Info<< "planarRefiner: total refinement events="
            << nRefinementIterations_ << ", total unrefinement passes="
            << nUnrefinementPasses_ << endl;
    }
    else
    {
        Info<< "planarRefiner: refinement iteration " << nRefinementIterations_
            << "/" << maxRefinementIterations_ << endl;
    }

    return changed;
}


void Foam::fvMeshTopoChangers::planarRefiner::topoChange
(
    const polyTopoChangeMap&
)
{
    // refinementLevel_ is a registered volScalarField and is mapped by
    // fvMesh::topoChange(). The reversible cutter is updated explicitly by
    // planarRefinementIterator/unrefineOnePass immediately after each map.
}


void Foam::fvMeshTopoChangers::planarRefiner::mapMesh(const polyMeshMap&)
{
    if (automaticUnrefinement_)
    {
        FatalErrorInFunction
            << "planarRefiner automaticUnrefinement v1.1.0 cannot map "
            << "reversible undo history across an arbitrary mesh map. Disable "
            << "automaticUnrefinement before mesh-to-mesh mapping."
            << exit(FatalError);
    }
}


void Foam::fvMeshTopoChangers::planarRefiner::distribute
(
    const polyDistributionMap&
)
{
    if (automaticUnrefinement_)
    {
        FatalErrorInFunction
            << "planarRefiner automaticUnrefinement v1.1.0 does not "
            << "redistribute undoableMeshCutter ancestry. Dynamic load "
            << "balancing/mesh redistribution must be disabled while "
            << "automaticUnrefinement is active." << exit(FatalError);
    }

    Info<< "planarRefiner: received native mesh-distribution callback" << endl;
}


bool Foam::fvMeshTopoChangers::planarRefiner::write(const bool doWrite) const
{
    if (!doWrite) return true;

    label nActiveReversibleSplits=0;
    if (automaticUnrefinement_ && reversibleCutter_.valid())
    {
        nActiveReversibleSplits=reversibleCutter_->getSplitFaces().size();
        reduce(nActiveReversibleSplits,sumOp<label>());
    }

    state_.set("formatVersion",2);
    state_.set("nRefinementIterations",nRefinementIterations_);
    state_.set("nUnrefinementPasses",nUnrefinementPasses_);
    state_.set("globalCells",mesh().globalData().nTotalCells());
    state_.set("geometry",geometry_);
    state_.set("automaticUnrefinement",bool(automaticUnrefinement_));
    state_.set("activeReversibleSplits",nActiveReversibleSplits > 0);

    state_.instance()=mesh().time().name();
    const bool ok=state_.regIOobject::write(doWrite);

    if (ok)
    {
        Info<< "planarRefiner: wrote restart state at time "
            << mesh().time().name() << " refinementEvents="
            << nRefinementIterations_ << " unrefinementPasses="
            << nUnrefinementPasses_;
        if (automaticUnrefinement_)
        {
            Info<< " activeReversibleSplits=" << nActiveReversibleSplits;
        }
        Info<< endl;
    }
    return ok;
}
