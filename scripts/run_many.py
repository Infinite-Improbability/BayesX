#! /usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path
from subprocess import run

# Parse command line arguments
parser = ArgumentParser(description="Run infile multiple times")
parser.add_argument(
    "infile",
    type=Path,
    help="Infile to run",
)
parser.add_argument("n_runs", type=int)
args = parser.parse_args()

infile_path: Path = args.infile

new_infile = []

out_path = "/dev/null"

with infile_path.open() as f:
    read_next = False
    for line in f:
        if read_next:
            out_path = Path(line.strip("'\n"))
            continue
        elif "#root" in line:
            read_next = True
            continue
        new_infile.append(line)

for i in range(1, args.n_runs + 1):
    iteration_path = Path(f"/tmp/{infile_path.stem}_it{i}.inp")
    iteration_infile = new_infile + ["#root\n", f"'{out_path}_it{i}'\n"]
    with iteration_path.open("x") as f:
        f.writelines(iteration_infile)

    run(
        [
            "mpiexec",
            "bin/BayesX",
            iteration_path,
        ],
        check=True,
    )
