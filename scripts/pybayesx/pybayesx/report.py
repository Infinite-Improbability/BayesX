from dataclasses import dataclass
from datetime import datetime
from itertools import zip_longest
from logging import getLogger
from pathlib import Path
from subprocess import CalledProcessError, run
from typing import Optional, Union

from .model import Model, Property, models
from .plot import plot

log = getLogger(__name__)


@dataclass
class Document:
    title: str
    priors: dict[str, tuple[str, str, str]]
    log_evidence: tuple[str, str]
    model: Model
    posteriors: dict[int, tuple[str, bool, str, str, str, str]]
    posterior_plot: Path
    parameters: dict[str, str]
    chains_path: Path
    comments: Optional[list[str]] = None

    def to_markdown(self, for_conversion_to=None):
        lines: list[str] = []

        # Special formatting for export types
        if for_conversion_to is not None:
            lines.append("---")
            lines.append(f"title: {self.title}")
            lines.append("---")

            lines.append("``` {=html}")
            lines.append("<style>body {max-width: unset} </style>")
            lines.append("```")

        # Title
        lines.append(f"# {self.title}")
        lines.append("\n")

        # Report generation time
        now = datetime.now().isoformat()
        lines.append(f"Report generated at {now}")
        lines.append("\n")

        lines.append("## Summary")
        free_priors = [
            p for p in self.model.required_priors if self.priors[p.value][0] != "Delta"
        ]
        lines.append(f"Model `{self.model.name}` used with {len(free_priors)} priors.")
        lines.append(
            f"Nested sampling global log-evidence is {' ± '.join(self.log_evidence)}."
        )
        lines.append(f"Chains saved to `{self.chains_path}`")
        lines.append("\n")

        # Comments
        if self.comments:
            lines.append("## Comments")
            lines += self.comments
            lines.append("\n")

        # Key parameters
        lines.append("## Key Parameters")
        key_params = ["nx", "ny", "xrayNch", "xrayEmin", "xrayEmax", "cluster_model"]
        lines += _make_md_table(
            ["Parameter", "Value"],
            key_params,
            [self.parameters[k] for k in key_params],
        )

        # Priors
        # TODO: Units
        lines.append(" ## Priors")
        # we want an ordered list so the columns stay ordered
        prior_keys = list(self.priors.keys())
        lines += _make_md_table(
            ["Prior", "Type", "Value 1", "Value 2"],
            [k for k in prior_keys],
            [self.priors[k][0] for k in prior_keys],
            [self.priors[k][1] for k in prior_keys],
            [self.priors[k][2] for k in prior_keys],
        )
        lines.append("For delta priors both values represent the function value.")
        lines.append(
            "For uniform and uniform-log priors value 1 is min and value 2 is max."
        )
        lines.append("For normal priors value 1 is mean and value 2 is mean.")
        lines.append("For log-normal priors value 1 is mean and value 2 is width.")
        lines.append("\n")

        # Posterior plot
        lines.append("## Posterior Plot")
        lines.append(f"![Plot of independent posteriors]({self.posterior_plot})")
        lines.append("\n")

        # Full input
        lines.append("## Full Parameters")
        parameter_keys = list(self.parameters.keys())
        parameter_keys.sort()
        lines += _make_md_table(
            ["Parameter", "Value"],
            parameter_keys,
            [self.parameters[k] for k in parameter_keys],
        )

        # Full posteriors
        lines.append("## All Posteriors")
        posterior_keys = list(self.posteriors.keys())

        if for_conversion_to == ".pdf" or for_conversion_to == ".tex":
            for k, v in self.posteriors.items():
                v = [f"${v[0]}$", v[1]] + [f"{float(n):.5e}" for n in v[2:]]
                self.posteriors[k] = tuple(v)

        lines += _make_md_table(
            [
                "Parameter",
                "Derived",
                "Mean",
                "Sigma",
                "Maximum Likelihood",
                "Maximum-A-Posteriori",
            ],
            [self.posteriors[k][0] for k in posterior_keys],
            [self.posteriors[k][1] for k in posterior_keys],
            [self.posteriors[k][2] for k in posterior_keys],
            [self.posteriors[k][3] for k in posterior_keys],
            [self.posteriors[k][4] for k in posterior_keys],
            [self.posteriors[k][5] for k in posterior_keys],
        )

        return lines


def make_report(
    infile_path: Union[Path, str], report_path: Optional[Path] = None, use_M500=False
):
    infile = _load_infile(Path(infile_path))

    chains_path = Path(infile["root"])

    if not chains_path.parent.is_dir():
        raise ValueError(f"Chains path {chains_path.parent} does not exist")

    if report_path is None:
        report_path = Path(str(chains_path) + "report.md")

    if report_path.exists():
        raise ValueError(f"Report path {report_path} already exists.")

    if "title" in infile:
        title = infile["title"]
    else:
        title = str(Path(*chains_path.parts[-2:]))

    comments = infile.get("comments", "")

    for m in models.values():
        if m.num == int(infile["cluster_model"]):
            model = m
            break
    else:
        model = Model(
            int(infile["cluster_model"]), [p for p in Property], "Unknown Model"
        )
    priors = {
        k.value: _parse_prior_values(infile[k.value]) for k in model.required_priors
    }

    log_evidence, posteriors = _load_stats_file(chains_path)

    # Pick posteriors to plot
    indep_posteriors = [k for k in posteriors if posteriors[k][1] is False]
    plot_parameters = [f"p{k:03}" for k in indep_posteriors]

    # Replace M200 with M500
    if use_M500:
        for i, p in enumerate(indep_posteriors):
            if "M_{T,200}" in posteriors[p][0]:
                m500_index = i
                for p2 in posteriors:
                    if "M_{T,500}" in posteriors[p2][0]:
                        m500_index = p2
                        break
                else:
                    break
                plot_parameters[i] = f"p{m500_index:03}"

    plot_file = Path(report_path.parent, report_path.stem + "_fig.svg")
    plot(
        chains_path,
        plot_parameters,
        display=False,
        save=True,
        plot_file=plot_file,
        chain_labels=title,
    )

    doc = Document(
        title=title,
        parameters=infile,
        priors=priors,  # type: ignore
        model=model,
        chains_path=chains_path,
        log_evidence=log_evidence,  # type: ignore
        posteriors=posteriors,  # type: ignore
        posterior_plot=plot_file.relative_to(report_path.parent),
        comments=[comments],
    )

    convert_to = None
    if report_path.suffix == "":
        report_path = report_path.with_suffix(".md")
    if report_path.suffix != ".md":
        convert_to = report_path.suffix
        report_path = report_path.with_suffix(".tmp.md")

    lines = doc.to_markdown(convert_to)

    with report_path.open("x") as f:
        f.writelines([line + "\n" for line in lines])
        log.info(f"Report saved to {report_path}")

    if convert_to:
        old_path = report_path
        report_path = report_path.with_stem(report_path.stem[:-4])
        report_path = report_path.with_suffix(convert_to)

        args = [
            "pandoc",
            old_path,
            "-o",
            report_path,
            "-f",
            "commonmark_x",
            "--standalone",
            "--table-of-contents",
        ]

        if convert_to in [".pdf", ".tex"]:
            args += [
                "-V",
                "papersize:a4",
                "-V",
                "geometry:margin=2cm",
            ]
        # landscape looks bad but I ran into pandoc errors trying to rotate the
        # oversized table

        log.info(f"Converting from markdown to {convert_to} with Pandoc")
        log.debug(f"Pandoc args: {args}")

        try:
            run(
                args,
                check=True,
            )
        except CalledProcessError as e:
            log.error("Error converting report")
            log.info(f"Removing {old_path}")
            old_path.unlink()
            raise e

        log.info(f"Report saved to {report_path}")
        log.info(f"Removing {old_path}")
        old_path.unlink()

    return report_path


def _load_infile(infile: Path) -> dict[str, str]:
    data = {}
    with infile.open() as f:
        reading = False
        key = None
        for line in f.readlines():
            if line[0] == "#":
                key = line.strip()[1:]  # strip whitespace and hash
                reading = True
                continue
            elif reading:
                data[key] = line.strip().strip("'")
                key = None
            reading = False

    return data


def _load_stats_file(chains_root: Path):
    # pathlib joinpath would insert a / so we treat it as a string instead
    path = Path(str(chains_root) + "stats.dat")

    evidence = None
    posteriors: dict[int, list] = {}

    posterior_names = _load_paramnames_file(chains_root)

    # Posteriors should be [name, derived, mean, sigma, max likelihood, MAP]

    with path.open() as stats:
        for line in stats.readlines():
            if "Log-Evidence" in line:
                evidence = [li.strip() for li in line.split(":")[1].split("+/-")]
                continue
            if line.strip() == "":
                continue
            elif "Parameters" in line:
                continue
            elif "Dim No" in line:
                continue
            sp = line.split()
            sp[0] = int(sp[0])  # type: ignore

            if sp[0] in posteriors:
                posteriors[sp[0]] += sp[1:]
            else:
                posteriors[sp[0]] = list(posterior_names[sp[0]]) + sp[1:]

    return (evidence, posteriors)


def _load_paramnames_file(chains_root: Path):
    # pathlib joinpath would insert a / so we treat it as a string instead
    path = Path(str(chains_root) + ".paramnames")

    # [name, is_derived]
    posteriors: dict[int, tuple[str, bool]] = {}

    with path.open() as pnames:
        for line in pnames.readlines():
            derived = False
            parts = [li.strip() for li in line.split()]
            key = parts[0].strip("p").strip("0")
            if key[-1] == "*":
                derived = True
                key = key[:-1]
            value = parts[1]
            if "/" in line:
                unit = line.split("/")[1].strip()
                unit = unit.replace("\\mathrm", "")
                unit = unit.replace("\\odot", "sun")
                # unit = unit.replace('{', '')
                # unit = unit.replace('}', '')
                value = f"{value} ({unit})"
            posteriors[int(key)] = (value, derived)

    return posteriors


def _make_md_table(headers: list[str], *args) -> list[str]:
    table: list[str] = []

    if len(args) < 1:
        raise ValueError("No columns specified")
    elif len(headers) > len(args):
        raise ValueError("More headers than columns")
    elif len(headers) < len(args):
        # Pad with empty headers
        headers += [""] * (len(args) - len(headers))

    column_widths = []
    for i, col in enumerate(args):
        column_widths.append(max([len(str(c)) for c in col + [headers[i - 1]]]))

    headers = [h.ljust(w) for h, w, in zip(headers, column_widths)]

    table.append(f"| {' | '.join(headers)} |")
    table.append(f"| {' | '.join(['-'*w for w in column_widths])} |")
    for cells in zip_longest(*args, fillvalue=""):
        cells = [str(cell).ljust(w) for cell, w in zip(cells, column_widths)]
        table.append(f"| {' | '.join(cells)} |")

    table.append("\n")

    return table


def _parse_prior_values(value: str):
    values = [v.strip() for v in value.split()]
    # TODO: Can we do better than hardcoding?
    prior_types = {
        "0": "Delta",
        "1": "Uniform",
        "2": "Uniform in log",
        "3": "Normal",
        "4": "Log-Normal",
    }
    values[0] = prior_types[values[0]]

    return values
