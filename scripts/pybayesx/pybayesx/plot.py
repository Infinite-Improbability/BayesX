from collections.abc import Iterable
from datetime import datetime
from logging import getLogger
from pathlib import Path
from typing import Optional, Sequence, Union

from getdist import MCSamples, loadMCSamples, plots

log = getLogger(__name__)
getLogger("matplotlib").setLevel("INFO")


def plot(
    chain_paths: Union[Path, str, Sequence[Path | str]],
    parameters: Iterable[str | tuple[str, float] | tuple[str, None]],
    display: bool = False,
    save: bool = True,
    plot_file: Optional[Path] = None,
    chain_labels: Union[str, Sequence[str]] = [],
):
    """Plot chains

    :param chain_paths: Path, or paths, to chains, including filename prefix,
     e.g. 'chains/subfolder/out' for files named 'outFILENAME'
    :type chain_paths: Path or list of paths, or string(s)
    :param parameters: Parameter ids to plot, e.g. [p001]
    :type parameters: Iterable[str  |  tuple[str, float]  |  tuple[str, None]]
    :param display: Show interactive plot, defaults to False
    :type display: bool
    :param save: Export plot to file, defaults to False
    :type save: bool
    :param plot_file: Path to export plot as. If None (default) file is saved in first
     chains folder as an SVG with a name based on the timestamp.
    :type plot_file: Optional[Path], optional
    :param chain_labels: Labels for the chains. Used in the legend, must be the same
     length as chain_paths, defaults to []
    :type chain_labels: Union[str, Sequence[str]], optional
    :raises Exception: Parameter name definition file missing.
    """
    # Massage input

    # This allows us to accept a single path without the user having to wrap it in a
    # sequence
    if isinstance(chain_paths, str):
        chain_paths = [Path(chain_paths)]
    elif isinstance(chain_paths, Path):
        chain_paths = [chain_paths]

    chain_paths = [Path(p) for p in chain_paths]

    log.debug(f"Plotting chains {chain_paths}")

    if isinstance(chain_labels, str):
        chain_labels = [chain_labels]

    # Set default output path
    if plot_file is None:
        plot_file = Path(
            f"{chain_paths[0]}_{datetime.now().strftime('%Y%m%d%H%M%S')}_tri.svg"
        )

    # Convert all parameters to tuples of the form (id_number, true_value) using None
    # if true value not provided
    parameters = [p if isinstance(p, tuple) else (p, None) for p in parameters]

    log.debug(f"Parmeters and values are: {parameters}")

    # Load in chains
    chains: list[MCSamples] = []
    for p in chain_paths:
        chains.append(loadMCSamples(str(p)))

    # Load parameter definitions for first chain
    first_chain_pars = chains[0].getParamNames()
    if first_chain_pars is None:
        raise Exception(
            f"MCSample {chains[0]} has no parameter names."
            "Is there a .paramnames file alongside the chains?"
        )

    # Select the parameters we are plotting
    plotted_pars: list[str] = []
    for p in parameters:
        plotted_pars.append(first_chain_pars.parWithName(p[0], error=True).label)

    # Internal names for the various parameter plots
    plotpar_names = ["plot" + str(i) for i in range(len(plotted_pars))]

    # Convert pairs of param id numbers and true values to dict
    markers = {
        p: v[1] for p, v in zip(plotted_pars, parameters, strict=True) if v is not None
    }

    log.debug(f"Markers are: {markers}")

    # Loop over all the chains
    for i, samps in enumerate(chains):
        parsi = samps.getParamNames()
        if parsi is None:
            raise Exception(
                f"MCSample {chains[0]} has no parameter names."
                "Is there a .paramnames file alongside the chains?"
            )

        # Search for parameters by name in current chain
        for p_label in plotted_pars:
            found = False

            for pj in parsi.names:
                # Found parameter
                if pj.label == p_label:
                    found = True

                    # Get internal name for parameter plot
                    p_name = plotpar_names[plotted_pars.index(p_label)]

                    # Handle special case parameters
                    if p_label[0] == "M" and "M_{\\odot}" in p_label:
                        # Adds derived parameter for mass to get around labelling issues
                        M = getattr(samps.getParams(), pj.name)
                        M1 = samps.ranges.getLower(pj.name)
                        M2 = samps.ranges.getUpper(pj.name)
                        p_label_new = p_label.replace(
                            "M_{\\odot}", "10^{14} M_{\\odot}"
                        )
                        samps.addDerived(M * 1e-14, name=p_name, label=p_label_new)
                        samps.updateBaseStatistics()
                        # Adding the derived parameter doesn't retain the range info
                        if M1 is None:
                            M1 = "N"
                        else:
                            M1 *= 1e-14
                        if M2 is None:
                            M2 = "N"
                        else:
                            M2 *= 1e-14
                        samps.setRanges({p_name: [M1, M2]})
                        if p_label in markers and i == 0:
                            markers[p_label] *= 1e-14
                    elif p_label[0] == "S" and "Jy" in p_label:
                        # Convert source fluxes to mJy
                        S = getattr(samps.getParams(), pj.name)
                        S1 = samps.ranges.getLower(pj.name)
                        S2 = samps.ranges.getUpper(pj.name)
                        samps.addDerived(
                            S * 1e3, name=p_name, label=p_label.replace("Jy", "mJy")
                        )
                        samps.updateBaseStatistics()
                        # Adding the derived parameter doesn't retain the range info
                        if S1 is None:
                            S1 = "N"
                        else:
                            S1 *= 1e3
                        if S2 is None:
                            S2 = "N"
                        else:
                            S2 *= 1e3
                        samps.setRanges({p_name: [S1, S2]})
                        if p_label in markers and i == 0:
                            markers[p_label] *= 1e3
                    else:
                        pj.name = p_name

            if not found:
                log.warn(f"Warning: {p_label} not found in {chain_paths[i]}")

        samps.setParamNames(parsi)  # This looks like it undoes all our changes

    # Setup for plots
    if chain_labels:
        legend = chain_labels
    else:
        legend = []
        for c in chains:
            legend.append(c.getName())

    plotter = plots.get_subplot_plotter(width_inch=8)
    plotter.settings.axes_fontsize = 8
    plotter.settings.alpha_filled_add = 0.4

    # Plot all the parameters
    plotter.triangle_plot(
        chains,
        plotpar_names,
        legend_labels=legend,
        filled=True,
        legend_loc="upper right",
    )

    # This should always be true after triangle_plot call
    # but pyright appreciates us confirming it. No harm in checking.
    assert plotter.subplots is not None

    # Add true values, if they exist
    for ipar, p in enumerate(plotted_pars):
        log.debug(f"Iterating over marker ({ipar}, {p})")
        if p in markers:
            # Add line for true value
            log.debug(f"Trying to add marker for {p}")
            ax = plotter.subplots[ipar, ipar]
            ax.axvline(x=markers[p], color="k")
            # Add star for true values in corrolation plots
            for jpar, p2 in enumerate(plotted_pars[ipar + 1 :]):
                if p2 in markers:
                    jplot = ipar + jpar + 1
                    ax = plotter.subplots[jplot, ipar]
                    ax.plot(markers[p], markers[p2], "*k")

    # Save as image
    if save:
        # TODO: Clobber protection
        plotter.export(str(plot_file))
        log.info(f"Saving fig to {plot_file}")

    if display:
        if plotter.fig is not None:
            log.info("Displaying fig")
            plotter.fig.show()
        else:
            log.error("Plot display failed. Figure does not exist.")
