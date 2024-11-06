import subprocess
import MDAnalysis as mda


def run_analysis():
    "./scripts/run_inference.py inference.output_prefix=example_outputs/design_partialdiffusion inference.input_pdb=input_pdbs/    2KL8.pdb 'contigmap.contigs=[79-79]' inference.num_designs=10 diffuser.partial_T=10"
    pass


if __name__ == "__main__":
    run_analysis()
