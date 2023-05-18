#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 14:51:09 2023

@author: Paolo Cozzi <bunop@libero.it>

Parse NeEstimator LD (tabular) results
"""

import re
import csv
import sys
import logging

from typing import Tuple
from pathlib import Path, PosixPath
from collections import namedtuple

logger = logging.getLogger(__name__)

loci_pattern = re.compile("Number of Loci = ([0-9]+)")
step_pattern = re.compile("individuals_([0-9]+)_step")

header = [
    'step',
    'n_loci',
    'population',
    'samp_size',
    'weighted_mean',
    'ind_alleles',
    'r2',
    'exp_r2',
    'ne',
    'parametric_ci_low',
    'parametric_ci_high',
    'jackknife_ci_low',
    'jackknife_ci_high',
    'eff_df'
]

NexLD = namedtuple("NexLD", header)


def skip_lines(handle, skip) -> Tuple[int, list]:
    logger.info(f"Skipping {skip} lines")

    # track skipped lines
    skipped = list()

    for i in range(skip):
        line = handle.readline().strip()
        position = handle.tell()

        logger.warning(f"Skipping: {line}")
        skipped.append(line)

    return position, skipped


def search_loci(handle):
    global loci_pattern

    match = re.search(loci_pattern, handle.readline())

    if match:
        return int(match.groups()[0])


def search_step(data_file: PosixPath):
    global step_pattern

    name = data_file.name
    match = re.search(step_pattern, name)

    if match:
        return int(match.groups()[0])


def process_datafile(data_file: PosixPath):
    # get step
    step = search_step(data_file)

    with open(data_file) as handle:
        # read number of loci
        position, skipped = skip_lines(handle, 4)
        n_loci = search_loci(handle)

        # skip until header
        position, skipped = skip_lines(handle, 9)

        # define a reader object
        reader = csv.reader(handle, delimiter='\t', lineterminator="\n")

        for line in reader:
            if not line:
                break

            line = [column.strip() for column in line]

            # some record are integers, some float
            for i in [1, 3, 11]:
                line[i] = int(line[i])

            for i in [2, 4, 5, 6, 7, 8, 9, 10]:
                line[i] = float(line[i])

            record = NexLD._make([step, n_loci] + line)
            yield record


if __name__ == "__main__":
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    data_files = Path('.').glob('*NexLD.txt')

    writer = csv.writer(sys.stdout, delimiter=",", lineterminator="\n")
    writer.writerow(header)

    for data_file in data_files:
        for record in process_datafile(data_file):
            writer.writerow(list(record))
