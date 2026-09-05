# Gate G2 status

**PASSED** — compact OF8 ↔ OF14 equivalence/profile release qualification passed within the laptop budget (900 s; elapsed 31 s).

The previously closed/passed 2400-cell Gate C3/A3 spatial-profile comparison is checksum-locked and passes the tighter G2 release envelope. The current OpenFOAM 14 trajectory was replayed to 2.7e-8 s and its mesh, shock position, primary fields, all 33 species fields, and sum(Y) closure were checked.

Fresh OF8 executable rerun: not requested. A fresh OF8 rebuild is not a mandatory pre-release gate because the accepted cross-version profile evidence is preserved and the production OF14 solver path has remained qualified; full-resolution cross-version studies remain deferred under POST_RELEASE_DEFERRED_QUALIFICATION.md.
