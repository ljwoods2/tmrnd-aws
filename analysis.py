import subprocess
import zarrtraj
import MDAnalysis as mda
from pathlib import Path

# conda activate mda && python analysis.py


def run_analysis():
    dir = Path("/root/inputs")
    filenames = [file for file in dir.iterdir() if file.is_file()]
    if len(filenames) != 1:
        raise ValueError
    file = filenames[0].resolve()

    u = mda.Universe(file, topology_format="PDB")

    n_res = len(u.residues)

    rf = [
        "docker",
        "run",
        "-i",
        "--rm",
        "--runtime=nvidia",
        "--gpus",
        "all",
        "-v",
        f"/root/models:/root/models",
        "-v",
        f"/root/inputs:/root/inputs",
        "-v",
        f"/root/outputs:/root/outputs",
        "rfdiffusion",
        f"inference.output_prefix=/root/outputs/motifscaffolding",
        f"inference.model_directory_path=/root/models",
        f"inference.input_pdb={file}",
        "inference.num_designs=1",
        "diffuser.T=15",
        f"contigmap.contigs=[{n_res}-{n_res}]",
    ]

    subprocess.run(rf, shell=False, check=True)

    new_prot = find_pdb_file("/root/outputs")

    if new_prot is None:
        raise ValueError

    new_prot = new_prot.resolve()

    commands = [
        [
            "gmx",
            "pdb2gmx",
            "-f",
            new_prot,
            "-o",
            "processed_protein.gro",
            "-p",
            "topol.top",
            "-ff",
            "amber03",
            "-water",
            "spc",
        ],
        [
            "gmx",
            "editconf",
            "-f",
            "processed_protein.gro",
            "-o",
            "protein_box.gro",
            "-c",
            "-d",
            "1.0",
            "-bt",
            "cubic",
        ],
        [
            "gmx",
            "solvate",
            "-cp",
            "protein_box.gro",
            "-o",
            "solvated.gro",
            "-p",
            "topol.top",
        ],
        [
            "gmx",
            "grompp",
            "-f",
            "em.mdp",
            "-c",
            "solvated.gro",
            "-p",
            "topol.top",
            "-o",
            "em.tpr",
        ],
        ["gmx", "mdrun", "-v", "-deffnm", "em"],
        [
            "gmx",
            "grompp",
            "-f",
            "npt.mdp",
            "-c",
            "em.gro",
            "-r",
            "em.gro",
            "-p",
            "topol.top",
            "-o",
            "npt.tpr",
        ],
        ["gmx", "mdrun", "-v", "-deffnm", "npt"],
        ["gmx", "editconf", "-f", "solvated.gro", "-o", "output.pdb"],
    ]

    # Execute each command in sequence
    for cmd in commands:
        result = subprocess.run(cmd, cwd="/root/gromacs", check=True)
    # /root/gromacs/solvated.gro

    command = [
        "aws",
        "s3",
        "cp",
        "/root/gromacs/output.pdb",
        f"s3://tmrnd-prototype/{file.name}_diffused.pdb",
    ]

    # Run the command without shell=True
    subprocess.run(command, check=True)

    u_final = mda.Universe(
        "/root/gromacs/solvated.gro", "/root/gromacs/npt.trr"
    )
    with mda.Writer(
        f"s3://tmrnd-prototype/{file.name}_diffused.zarrmd",
        n_atoms=u_final.trajectory.n_atoms,
    ) as W:
        for ts in u_final.trajectory:
            W.write(u_final.atoms)


def find_pdb_file(directory):
    # Iterate through all files with .pdb extension recursively
    for file in Path(directory).rglob("*.pdb"):
        return file  # Return the first PDB file found

    return None  # If no PDB file is found


if __name__ == "__main__":
    run_analysis()
