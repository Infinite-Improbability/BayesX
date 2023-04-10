from argparse import ArgumentParser, Namespace
from pathlib import Path

from getdist import loadMCSamples

from .plot import plot


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
        "-i", "--interactive", type=bool, default=True, help="Show plot"
    )
    plot_parser.add_argument(
        "-l", "--label", type=str, nargs="+", default=[], help="Labels for chains"
    )
    plot_parser.set_defaults(func=interactive_plot)

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
        display=args.interactive,
        plot_file=args.output,
        chain_labels=args.label,
    )
