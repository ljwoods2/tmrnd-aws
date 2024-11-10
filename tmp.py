import zarrtraj
import MDAnalysis as mda
import fsspec

with fsspec.open(
    "s3://tmrnd-prototype/sample_job_diffused.pdb", "r", anon=True
) as top:
    storage_options = {"anon": True}
    u = mda.Universe(
        top,
        "s3://tmrnd-prototype/sample_job_diffused.zarrmd",
        storage_options=storage_options,
        topology_format="PDB",
    )
    print(u)

import zarrtraj
import MDAnalysis as mda
import fsspec
import MDAnalysis.analysis.rms

with fsspec.open(
    "s3://tmrnd-prototype/sample_job_diffused.pdb", "r", anon=True
) as top:
    storage_options = {"anon": True}
    u = mda.Universe(
        top,
        "s3://tmrnd-prototype/sample_job_diffused.zarrmd",
        storage_options={"anon": True},
        topology_format="PDB",
    )
    R = MDAnalysis.analysis.rms.RMSD(u, u, select="backbone")
    R.run()
    results = R.rmsd
    print(results)
