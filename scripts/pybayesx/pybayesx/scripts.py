import logging
from argparse import ArgumentParser, Namespace
from pathlib import Path
from sys import stderr, stdout

from getdist import loadMCSamples

from .analysis import demo
from .mask import mask
from .plot import plot
from .report import make_report

log = logging.getLogger(__name__)


def main():
    parser = ArgumentParser(description="Basic interaction with BayesX")
    subparsers = parser.add_subparsers(required=True, help="sub-command")

    # Plot subparser
    plot_parser = subparsers.add_parser(
        "plot", description="Generate plots from chains"
    )
    plot_parser.add_argument(
        "chain_paths", type=Path, help="Paths to chains", nargs="+"
    )
    plot_parser.add_argument(
        "-p" "--parameter",
        type=str,
        nargs="+",
        help="Parameter ids to plot",
        default=[],
    )
    plot_parser.add_argument(
        "-o", "--output", type=Path, default=None, help="Path to save plot as"
    )
    plot_parser.add_argument(
        "-ni", "--no-interactive", action="store_false", help="Do not display plot"
    )
    plot_parser.add_argument(
        "-l", "--label", type=str, nargs="+", default=[], help="Labels for chains"
    )
    plot_parser.set_defaults(func=interactive_plot)

    # Mask subparser
    mask_parser = subparsers.add_parser(
        "mask",
        description="Generate mask for given files of ellipses generated by CIAO",
    )
    bounds = mask_parser.add_argument_group("bounds", "Coordinates of data bounds.")
    bounds.add_argument("xMin", type=float, help="Minimum x coordinate in data.")
    bounds.add_argument("yMin", type=float, help="Minimum y coordinate in data.")
    bounds.add_argument("xMax", type=float, help="Maximum x coordinate in data.")
    bounds.add_argument("yMax", type=float, help="Maximum x coordinate in data.")
    mask_parser.add_argument(
        "maskFiles", nargs="+", help="List of paths to files of ellipses to mask"
    )
    mask_parser.set_defaults(func=_mask)

    # Demo parser
    demo_parser = subparsers.add_parser("demo", description="Run demo analysis.")
    demo_parser.set_defaults(func=_demo)

    # Report parser
    report_parser = subparsers.add_parser(
        "report", description="Write a report on a run."
    )
    report_parser.add_argument("infile", type=Path, help="Path to infile")
    report_parser.add_argument(
        "-o", "--output", type=Path, default=None, help="Path to save report as"
    )
    report_parser.add_argument(
        "--m500", action="store_true", help="Use M500 instead of M200 on plot"
    )
    report_parser.set_defaults(func=_report)

    # Write logging to shell
    rlog = logging.getLogger()
    rlog.setLevel(logging.INFO)  # log everything
    infoHandler = logging.StreamHandler(stream=stdout)
    infoHandler.setLevel(logging.INFO)  # write info to stdout
    rlog.addHandler(infoHandler)
    errorHandler = logging.StreamHandler(stream=stderr)
    errorHandler.setLevel(logging.WARNING)  # warnings or errors go to stderr
    rlog.addHandler(errorHandler)

    # Call appropriate function for subcommand
    args = parser.parse_args()
    args.func(args)


def interactive_plot(args: Namespace):
    if not args.p__parameter:
        if len(args.chain_paths) > 1:
            print(
                "Select parameters based on first chain entry; "
                "will look for matching parameter names in subsequent entries"
            )

        ch = loadMCSamples(str(args.chain_paths[0]))
        pars = ch.getParamNames()

        if pars is None:
            raise Exception(
                f"MCSample {ch} has no parameter names."
                "Is there a .paramnames file alongside the chains?"
            )

        print("Parameters are: ")
        print(pars)
        plotpars = []
        while True:
            p = input("Parameter to plot or q to quit: ").strip()
            if p == "q":
                break
            if pars.parWithName(p) is None:
                print("Parameter not found")
            else:
                plotpars.append(p)

        args.p__parameter = plotpars

    if not args.output:
        default_output = Path("plot_tri.svg")
        ouput = input(f"Path to output file ({default_output}): ")
        args.output = Path(ouput) if ouput != "" else default_output

    # TODO: Could add an option to suppress export

    plot(
        args.chain_paths,
        args.p__parameter,
        display=args.no_interactive,
        plot_file=args.output,
        chain_labels=args.label,
    )


def _mask(args):
    mask(
        args.xMin,
        args.xMax,
        args.yMin,
        args.yMax,
        args.maskFiles,
    )


def _demo(args):
    demo()


def _report(args):
    rp = make_report(args.infile, args.output, args.m500)
