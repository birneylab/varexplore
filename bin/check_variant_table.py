#!/usr/bin/env python


"""
Provide a command line tool to validate and transform tabular samplesheets.
"""


import argparse
import csv
import logging
import sys
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a
            previously validated and transformed row.
            The order of rows is maintained.
    """

    def __init__(
        self,
        var_col="var_id",
        first_col="chr",
        second_col="region_start",
        third_col="region_end",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            var_col (str): The name of the column that contains the variant ID
                name (default "var_id").
            first_col (str): The name of the column that contains the chromosome ID
                path (default "chr").
            second_col (str): The name of the column that contains the starting position of the region of interest (default "region_start").
            third_col (str): The name of the column that contains the ending position of the region of interest (default "region_end").
        """
        super().__init__(**kwargs)
        self._var_col = var_col
        self._first_col = first_col
        self._second_col = second_col
        self._third_col = third_col
        self._seen = {
            self._var_col: set(),
            self._first_col: set(),
            self._second_col: set(),
            self._third_col: set(),
        }
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row.

        Args:
            row (dict): A mapping from column headers (keys) to elements of
                that row (values).

        """
        self._validate_var(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_third(row)
        self._seen[self._var_col].add(row[self._var_col])
        self._seen[self._first_col].add(row[self._first_col])
        self._seen[self._second_col].add(row[self._second_col])
        self._seen[self._third_col].add(row[self._third_col])
        self.modified.append(row)

    def _validate_var(self, row):
        """
        Assert that the sample name exists and convert spaces to underscores.
        """

        if len(row[self._var_col]) <= 0:
            raise AssertionError("Variant ID is required.")
        # Sanitize samples slightly.
        row[self._var_col] = row[self._var_col].replace(" ", "_")

    def _validate_first(self, row):
        """Check that the chromosome ID column is present"""

        if len(row[self._first_col]) <= 0:
            raise AssertionError("The chromosome ID is required.")

    def _validate_second(self, row):
        """Check that the region start column is present"""

        if len(row[self._second_col]) <= 0:
            raise AssertionError("The region start column is required.")

        if not (str(row[self._second_col]).isnumeric() and int(row[self._second_col]) > 0):
            raise AssertionError(
                "The region start entries must be integers greather than 0",
                f"Error detected for var_id: {row[self._var_col]}, region_start: {row[self._second_col]}"
            )

    def _validate_third(self, row):
        """Check that the region start column is present"""

        if len(row[self._third_col]) <= 0:
            raise AssertionError("The region end column is requries.")

        if not (str(row[self._third_col]).isnumeric() and int(row[self._third_col]) > int(row[self._second_col])):
            raise AssertionError(
                "The region end entries must be integers greather than region_start",
                f"Error detected for var_id: {row[self._var_col]}, region_start: {row[self._second_col]}, region_end: {row[self._third_col]}"
            )

    def validate_unique_var(self):
        """
        Assert that the sample namess are all unique.
        """

        if len(self._seen[self._var_col]) != len(self.modified):
            raise AssertionError("Variant IDs must be unique.")

def read_head(handle, num_lines=10):
    """
    Read the specified number of lines from the current position in the file.
    """
    lines = []

    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)

    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read
        position is expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)

    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core
    pipelines.

    Validate the general shape of the table, expected columns, and each row.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can
            be either CSV, TSV, or any other format automatically recognized by
            ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed
            samplesheet should be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following
        structure

            sample,cram,crai
            SAMPLE,SAMPLE.cram,SAMPLE.cram.crai
            SAMPLE,SAMPLE.cram,SAMPLE.cram.crai
    """
    required_columns = {"var_id", "chr", "region_start", "region_end"}
    optional_columns = set()
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on
    # `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.

        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(
                "The sample sheet **must** contain these column headers:"
                + f" {req_cols}."
            )
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()

        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} on line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_var()
    header = list(required_columns | optional_columns)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on
    # `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()

        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog=(
            "Example: python check_samplesheet.py samplesheet.csv"
            + " samplesheet.valid.csv"
        ),
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )

    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(
        level=args.log_level, format="[%(levelname)s] %(message)s"
    )

    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
