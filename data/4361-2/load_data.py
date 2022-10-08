import pybayesx.pybayesx as pb
from datetime import datetime

userlabel = input("Enter label for processed data: ")
b = pb.BayesX(
    f"4361-{userlabel}-{datetime.now().strftime('%Y%m%d%H%M')}",
    chain_root="chains/",
)
b.make_config()
b.set_priors()
b.load_arf("data/4361-2/4361_er.arf")
b.load_rmf("data/4361-2/4361_er.rmf")
b.load_events_from_fits("data/4361-2/4361_er_evts.fits", mode="evts")
b.load_events_from_fits("data/4361-2/4361_er_bg.fits", mode="bg")
b.load_mask("data/4361-2/mask.npz")

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
