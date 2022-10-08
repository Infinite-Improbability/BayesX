import pybayesx.pybayesx as pb
from datetime import datetime

userlabel = input("Enter label for processed data: ")
b = pb.BayesX(
    f"TNG3001s91h0-{userlabel}-{datetime.now().strftime('%Y%m%d%H%M')}",
    chain_root="chains/",
)
b.make_config()
b.set_priors()
b.load_arf("data/tngdata/acisi_aimpt_cy0.arf")
b.load_rmf("data/tngdata/acisi_aimpt_cy0.rmf")
b.load_events_from_fits("data/tngdata/tng_obs_evt.fits", mode="evts")
b.load_events_from_fits("data/tngdata/tng_bg.fits", mode="bg")
while True:
    n_bins = int(input("Select number of bins: "))
    cellsize = int(input("Select cellsize: "))
    b.bin_and_export(n_bins, cellsize)
    y = input("Enter y to accept: ")
    if y == "y":
        break

b.export_rmf()
b.export_arf()
b.export_config()
