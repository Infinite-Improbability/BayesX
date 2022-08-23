#! /usr/bin/env python3

"""
Automates internal consistency checking: whether BayesX can retrieve the parameters it used
to generate a prediction x-ray count.

# Usage
Run this script from the BayesX root folder.
It expects the BayesX executable to be located at bin/BayesX and compiled for parallel execution
with mpiexec.

# Process
1. Generate predicted x-ray counts using BayesX for fixed priors.
   Fixed priors are selected from a uniform distribution between prior.min and prior.max,
   which are specified in the prior definition below.
2. Apply poisson noise to predicted counts.
3. Free noisy prediction back into BayesX with free priors and the same model.
4. Generate plot of posteriors corresponding to free priors for easy analysis.
   True values are marked, if in plot range.

# Notes
* Chains are output to a subfolder the main chains directory with the datetime in the name.
  This is to keep things tidy and prevent successive runs from conflicting. This folder also
  includes the generated .inp files so as to record run configuration.
* Defaults to using all avaliable cores.
"""

import numpy as np
from subprocess import run
from os import cpu_count, mkdir
from datetime import datetime

rng = np.random.default_rng()


class Path:
    """
    Path string.

    Created because BayesX export wants paths in single quotes but
    in python we just want the strings. By using a custom __repr__ method
    single quotes are automatically applied on export without needing to check
    if it is a path
    """

    def __init__(self, path) -> None:
        self.path = path

    def __repr__(self) -> str:
        return f"'{self.path}'"


class Prior:
    """
    Define prior, using types defined in BayesX's readme.

    True value is selected from a uniform distribtion.
    TODO: Add support for other distributions.
    """

    def __init__(self, type_: int, min_: int, max_: int):
        self.type = type_
        self.min = min_
        self.max = max_
        if self.type == 0:
            assert self.min == self.max
            self.value = self.min
        else:
            self.value = rng.uniform(self.min, self.max)

    def free(self) -> str:
        return f"{self.type}    {self.min}  {self.max}"

    def fixed(self) -> str:
        return f"0    {self.value}  {self.value}"


# Get current time and use it to uniquely identify this run.
now = datetime.now().strftime("%Y%m%d%H%M%S")
chain_path = f"chains/auto-{now}"
mkdir(chain_path)

# All parameters will be defined in this dict
# Except priors, which have their own dict.
params = {}

# The root for MEKAL data files:
params["filion"] = Path("data/MEKAL/mekal1.dat")
params["filrec"] = Path("data/MEKAL/mekal2.dat")
params["filkar"] = Path("data/MEKAL/mekal3.dat")
params["filcon"] = Path("data/MEKAL/mekal4.dat")
params["fillin"] = Path("data/MEKAL/mekal5.dat")
params["filcaf"] = Path("data/MEKAL/mekal6.dat")

# The telescope
params["XrayTelescope"] = "CHA"  # The first three letters of the X-ray telescope name
params["filARF"] = Path(
    "data/simtestdatafiles/ARF_32by1.txt"
)  # Root for ARF telescope file in .txt format.
params["filRMF"] = Path(
    "data/simtestdatafiles/RMF_32by32.txt"
)  # Root for RMF telescope file in .txt format.

# Misc params
params["n"] = 64  # Number of steps for discretising r
params["nx"] = 64  # Number of pixels in x direction
params["ny"] = 64  # Number of pixels in y direction
params["xrayNbin"] = 32  # Number of energy bins
params["xrayNch"] = 32  # Number of energy bins
params["Aeffave"] = 250  # Average effective area of the telescope in cm^{2}
params["xraycell"] = 0.492
params["xrayEmin"] = 0.7
params["xrayEmax"] = 7.0
params["sexpotime"] = 120e3
params["bexpotime"] = 120e3
params["NHcol"] = 2.2e20
params["xrayBG_model"] = 8.4e-6
params["nlive"] = 100
params["eff"] = 0.8
params["tol"] = 0.5
params["seed"] = -1

# Currently only supports cluster model 1
params["cluster_model"] = 1

# Priors
# Note that the automatic free prior detection at the end assumes
# this specific set of prior types
priors = {}
priors["x_prior"] = Prior(0, 1, 1)
priors["y_prior"] = Prior(0, 1, 1)
priors["m200_prior"] = Prior(1, 4e14, 6e14)
priors["fgas200_prior"] = Prior(1, 0.1, 0.3)
priors["a_GNFW_prior"] = Prior(0, 1.062, 1.062)
priors["b_GNFW_prior"] = Prior(0, 5.4807, 5.4807)
priors["c_GNFW_prior"] = Prior(0, 0.3292, 0.3292)
priors["c500_GNFW_prior"] = Prior(0, 1.156, 1.156)
priors["alpha_model2_prior"] = Prior(1, 0.05, 5)
priors["gamma0_poly_prior"] = Prior(1, 0.1, 0.4)
priors["gammaR_poly_prior"] = Prior(1, 0.1, 0.4)
priors["t0_poly_prior"] = Prior(1, 1, 5)
priors["z_Prior"] = Prior(0, 0.0231, 0.0231)


# Start configuring for fixed priors
fixed_path = Path(f"{chain_path}/infile-auto-fixed.inp")
params["root"] = Path(f"{chain_path}/fixed_")
params["nCdims"] = "0"

# For the event files we'll use a blank background and mask
# For the fixed case we don't care about the event file, so we'll also use the blank file
blank = np.zeros(
    (params["nx"] * params["ny"] * params["xrayNch"], 1)
)  # BAYES-X reads the counts in a 1-D array of dimension of (nx times ny times xrayNch).
blank_path = Path(f'data/blank{params["nx"]}x{params["ny"]}x{params["xrayNch"]}.txt')
np.savetxt(blank_path.path, blank, fmt="%.1d")
params["filBG"] = blank_path  # The root for background file
params["filevent"] = blank_path  # The root for the event file
params["filmask"] = blank_path  # The root for the mask file.

# Export config to file
with open(fixed_path.path, "w") as f:
    for key, value in params.items():
        f.writelines(f"#{key}\n")
        f.write(str(value) + "\n")
        f.write("\n")
    for key, value in priors.items():
        f.write(f"#{key}\n")
        f.write(str(value.fixed()) + "\n")
        f.write("\n")

# And run for fixed case
time = datetime.now()
run(["mpiexec", "-n", str(cpu_count()), "bin/BayesX", fixed_path.path])
fixed_time = datetime.now() - time
print(f"\nFixed run complete in {fixed_time.total_seconds()} seconds\n\n")


# Load in results and apply poisson noise
fixed_data = np.loadtxt(params["root"].path + "generated-data.txt")
noisy_path = Path(params["root"].path + "generated-data-noisy.txt")
np.savetxt(noisy_path.path, rng.poisson(fixed_data), fmt="%.1d")

# Configure for free priors
free_path = Path(f"{chain_path}/infile-auto-free.inp")
params["root"] = Path(f"{chain_path}/free_")
params["nCdims"] = "2"
params["filevent"] = noisy_path

# And export to file
with open(free_path.path, "w") as f:
    for key, value in params.items():
        f.write(f"#{key}\n")
        f.write(str(value) + "\n")
    for key, value in priors.items():
        f.write(f"#{key}\n")
        f.write(str(value.free()) + "\n")
        f.write("\n")

# Run with free priors
time = datetime.now()
run(["mpiexec", "-n", str(cpu_count()), "bin/BayesX", free_path.path])
free_time = datetime.now() - time
print(f"\nFree run complete in {free_time.total_seconds()} seconds\n\n")

# Generate plot, with automatic detection of free priors
# TODO: Make this more robust
# Right now this will break, possibly silently, if which priors are fixed changes
plot_priors = []
true_priors = ["-t"]
p_count = 1
for p in priors.values():
    if p.type != 0:
        plot_priors.append(f"p{str(p_count).zfill(3)}")
        true_priors.append(str(p.value))
        p_count += 1
    # Break if we've addedd all the relevant params
    if params["cluster_model"] == 1 and p_count > 2:
        break
    elif params["cluster_model"] == 2 and p_count > 3:
        break
    elif params["cluster_model"] == 3 and p_count > 4:
        del plot_priors[1:2] # remove fgas200 and alpha_model2_prior
        del true_priors[1:2]
        break
    
run(
    ["python3", "scripts/auto_plot_tri.py", params["root"].path]
    + plot_priors
    + true_priors
)

# Output any additional data
with open(f"{chain_path}/autorun-report.txt", "x") as f:
    f.write(f"Fixed runtime: {fixed_time.total_seconds()} seconds\n")
    f.write(f"Free runtime: {free_time.total_seconds()} seconds")
    # TODO: Include priors, posterior information, etc
