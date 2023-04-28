"""Helper module for configuring parameters. Adapted from Kayran Schmidt."""
from __future__ import annotations

import argparse
import os
import re
from typing import Any, Sequence, TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import _FormatterClass


class Parameters(argparse.Namespace):
    """Namespace with regular expression matching for parameters.hpp file."""

    REGEX = {
        "kModel": re.compile(r"(#include \"models\/).+(\.hpp\")$", flags=re.M),
        "kR": re.compile(r"(Index kR = ).*(;)", flags=re.M),
        "kM1": re.compile(r"(Index kM1 = ).*(;)", flags=re.M),
        "kM2": re.compile(r"(Index kM2 = ).*(;)", flags=re.M),
        "kN1": re.compile(r"(Index kN1\[kM1\] {).*(};)", flags=re.M),
        "kN2": re.compile(r"(Index kN2\[kM2\] {).*(};)", flags=re.M),
        "kBinsize1": re.compile(r"(Index kBinsize1\[kM1\] {).*(};)", flags=re.M),
        "kBinsize2": re.compile(r"(Index kBinsize2\[kM2\] {).*(};)", flags=re.M),
        "kLiml1": re.compile(r"(double kLiml1\[kM1\] {).*(};)", flags=re.M),
        "kLiml2": re.compile(r"(double kLiml2\[kM2\] {).*(};)", flags=re.M),
        "kTstar": re.compile(r"(double kTstar = ).*(;)", flags=re.M),
        "kTau": re.compile(r"(double kTau = ).*(;)", flags=re.M),
        "kSnapshot": re.compile(r"(Index kSnapshot = ).*(;)", flags=re.M),
        "kSecondOrder": re.compile(r"(bool kSecondOrder = ).*(;)", flags=re.M),
        "kNSubsteps": re.compile(r"(Index kNSubsteps = ).*(;)", flags=re.M),
        "kFilename": re.compile(r"(char kFilename\[\] = \").*(\";)", flags=re.M)
    }

    def configure(self):
        """Read parameters.hpp and replace specified settings.

        Take parameters.hpp file located at path set by -o / --output.
        """
        read_path = write_path = os.path.join(self.OUTPUT, "parameters.hpp")

        with open(read_path, "r", encoding="utf8") as f:
            parameters_str = f.read()

        params_to_write = set(
            filter(lambda p: p[0] in self.REGEX, vars(self).items()))
        for param, val in params_to_write:
            # Convert Python to C++ syntax
            if isinstance(val, tuple):
                valstr = str(val)[1:-1]
            elif isinstance(val, bool):
                valstr = str(val).lower()
            else:
                valstr = str(val)
            parameters_str = re.sub(
                self.REGEX[param],
                r"\g<1>" + valstr + r"\g<2>",
                parameters_str
            )

        print(
            "\nThese parameters were set but will not be written to `parameters.hpp`:\n\t" +
            "\n\t".join(map(str, set(vars(self).items()) - params_to_write))
        )

        with open(write_path, "w", encoding="utf8") as fout:
            fout.write(parameters_str)


class ParametersParser(argparse.ArgumentParser):
    """Customized ArgumentParser for use with parameters.hpp file."""

    def __init__(
        self,
        prog: str | None = None,
        usage: str | None = None,
        description: str | None = None,
        epilog: str | None = None,
        parents: Sequence[argparse.ArgumentParser] = [],
        formatter_class: _FormatterClass = argparse.HelpFormatter,
        prefix_chars: str = '-',
        fromfile_prefix_chars: str | None = None,
        argument_default: Any = argparse.SUPPRESS,
        conflict_handler: str = 'error',
        add_help: bool = True,
        allow_abbrev: bool = True
    ) -> None:

        description += (
            "\nOptional arguments will only be written to file if they are set."
            " Unspecified optional arguments are taken from the existing file."
        )

        super().__init__(
            prog,
            usage,
            description,
            epilog,
            parents,
            formatter_class,
            prefix_chars,
            fromfile_prefix_chars,
            argument_default,
            conflict_handler,
            add_help,
            allow_abbrev
        )

        self.add_argument(
            '-o',
            '--output',
            default="../include",
            type=str,
            metavar="path/to/input",
            help="Path to folder for saving generated files. `parameters.hpp` has to be located in `project_root/include`, therefore the default value is set to `../include`. This has to be adjusted if the script is located in a directory other than `project_root/input`.",
            dest="OUTPUT")

    def parse_args(self, args=None, namespace: Parameters = None) -> Parameters:
        if namespace is None:
            namespace = Parameters()

        return super().parse_args(args, namespace)

    def add_model(self):
        self.add_argument(
            '--model',
            type=str,
            metavar="kModel",
            dest="kModel",
            help="Name of the model file (without suffix `.hpp`). Model files have to be located in `include/models/`."
        )

    def add_rank(self):
        self.add_argument(
            '--rank',
            type=int,
            metavar="kR",
            dest="kR",
            help="Rank of the singular value decomposition."
        )

    def add_dimension1(self):
        self.add_argument(
            '--dim1',
            type=int,
            metavar="kM1",
            dest="kM1",
            help="Number of species in partition 1."
        )

    def add_dimension2(self):
        self.add_argument(
            '--dim2',
            type=int,
            metavar="kM2",
            dest="kM2",
            help="Number of species in partition 2."
        )

    def add_gridpoints1(self):
        self.add_argument(
            '--n1',
            nargs="+",
            type=int,
            metavar="kN1",
            dest="kN1",
            help="Number of gridpoints (bins) for partition 1."
        )

    def add_gridpoints2(self):
        self.add_argument(
            '--n2',
            nargs="+",
            type=int,
            metavar="kN2",
            dest="kN2",
            help="Number of gridpoints (bins) for partition 2."
        )

    def add_binsize1(self):
        self.add_argument(
            '--bin1',
            nargs="+",
            type=int,
            metavar="kBinsize1",
            dest="kBinsize1",
            help="Binsize for population numbers in partition 1."
        )

    def add_binsize2(self):
        self.add_argument(
            '--bin2',
            nargs="+",
            type=int,
            metavar="kBinsize2",
            dest="kBinsize2",
            help="Binsize for population numbers in partition 2."
        )

    def add_liml1(self):
        self.add_argument(
            '--liml1',
            nargs="+",
            type=int,
            metavar="kLiml1",
            dest="kLiml1",
            help="Lower limits for population numbers in partition 1."
        )

    def add_liml2(self):
        self.add_argument(
            '--liml2',
            nargs="+",
            type=int,
            metavar="kLiml2",
            dest="kLiml2",
            help="Lower limits for population numbers in partition 2."
        )

    def add_tstar(self):
        self.add_argument(
            '--tstar',
            type=float,
            metavar="kTstar",
            dest="kTstar",
            help="Final time for integration."
        )

    def add_tau(self):
        self.add_argument(
            '--tau',
            type=float,
            metavar="kTau",
            dest="kTau",
            help="Time step size for integration."
        )

    def add_snapshot(self):
        self.add_argument(
            '--snapshot',
            type=int,
            metavar="kSnapshot",
            dest="kSnapshot",
            help="Number of time steps between a snapshot."
        )

    def add_secondorder(self):
        self.add_argument(
            '--so',
            type=bool,
            metavar="kSecondOrder",
            dest="kSecondOrder",
            help="Flag for activating the second order method."
        )

    def add_substeps(self):
        self.add_argument(
            '--substeps',
            type=int,
            metavar="kNSubsteps",
            dest="kNSubsteps",
            help="Number of substeps for second order method."
        )

    def add_filename(self):
        self.add_argument(
            '--fname',
            type=str,
            metavar="kFilename",
            dest="kFilename",
            help="Name of the output folder (input string without quotes). Folder will be located in `output`."
        )