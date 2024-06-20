#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 16:05:41 2024

@author: jkoubele
"""

import os

os.environ[
    'OPENBLAS_NUM_THREADS'] = '1'  # solves weird error when importing numpy (and consequently e.g. pandas, biopython etc.) on cluster

import argparse
import json
from pathlib import Path
from typing import Optional

import pysam
from pybedtools import BedTool, Interval
from interval import interval as py_interval


def extract_id_of_invalid_reads(bamfile_input_path: Path) -> set[str]:
    bamfile_input = pysam.AlignmentFile(bamfile_input_path, "rb")
    invalid_ids: set[str] = set()
    for i, read in enumerate(bamfile_input):
        if i % 100_000 == 0:
            print(f"Finding invalid reads: {i}", flush=True)
        if read.is_secondary:
            invalid_ids.add(read.query_name)
    bamfile_input.close()
    return invalid_ids


def extract_and_save_unique_pairs(bamfile_input_path: Path,
                                  write_bam_files=True,
                                  output_folder: Optional[Path] = None,
                                  output_bam_file_forward_path: Optional[Path] = None,
                                  output_bam_file_reverse_path: Optional[Path] = None,
                                  output_bed_file_forward_path: Optional[Path] = None,
                                  output_bed_file_reverse_path: Optional[Path] = None,
                                  output_json_read_count_file: Optional[Path] = None) -> None:
    output_folder = bamfile_input_path.parent if output_folder is None else output_folder
    output_bam_file_forward_path = output_folder / 'unique_forward.bam' \
        if output_bam_file_forward_path is None else output_bam_file_forward_path
    output_bam_file_reverse_path = output_folder / 'unique_reverse.bam' \
        if output_bam_file_reverse_path is None else output_bam_file_reverse_path
    output_bed_file_forward_path = output_folder / 'unique_forward.bed.gz' \
        if output_bed_file_forward_path is None else output_bed_file_forward_path
    output_bed_file_reverse_path = output_folder / 'unique_reverse.bed.gz' \
        if output_bed_file_reverse_path is None else output_bed_file_reverse_path
    output_json_read_count_file = output_folder / 'unique_read_count.json' \
        if output_json_read_count_file is None else output_json_read_count_file

    invalid_ids = extract_id_of_invalid_reads(bamfile_input_path)

    reads_1: dict[str, pysam.AlignedSegment] = {}
    reads_2: dict[str, pysam.AlignedSegment] = {}

    intervals_forward: list[Interval] = []
    intervals_reverse: list[Interval] = []

    bamfile_input = pysam.AlignmentFile(bamfile_input_path, "rb")

    valid_reads_forward = 0
    valid_reads_reverse = 0
    for i, read in enumerate(bamfile_input):
        if i % 100_000 == 0:
            print(f"Computing covered intervals: {i}", flush=True)
        if read.query_name in invalid_ids:
            continue
        if read.is_read1 and read.query_name not in reads_2:
            reads_1[read.query_name] = read
            continue

        elif read.is_read2 and read.query_name not in reads_1:
            reads_2[read.query_name] = read
            continue

        if read.is_read1 and read.query_name in reads_2:
            read_1 = read
            read_2 = reads_2.pop(read.query_name)

        elif read.is_read2 and read.query_name in reads_1:
            read_1 = reads_1.pop(read.query_name)
            read_2 = read
        else:
            assert False
        if read_1.reference_name != read_2.reference_name:
            invalid_ids.add(read_1.query_name)
            continue

        interval_union = py_interval(*(read_1.get_blocks() + read_2.get_blocks()))
        interval_union = sorted(list(interval_union), key=lambda x: x[0])
        if read_1.is_forward:
            intervals_forward.extend([Interval(chrom=read_1.reference_name,
                                               start=int(x[0]),
                                               end=int(x[1]),
                                               strand='+') for x in interval_union])
            valid_reads_forward += 2
        elif read_1.is_reverse:
            intervals_reverse.extend([Interval(chrom=read_1.reference_name,
                                               start=int(x[0]),
                                               end=int(x[1]),
                                               strand='-') for x in interval_union])
            valid_reads_reverse += 2
    bamfile_input.close()

    print("Saving intervals to .bed files", flush=True)
    BedTool(intervals_forward).sort().saveas(output_bed_file_forward_path, compressed=True)
    BedTool(intervals_reverse).sort().saveas(output_bed_file_reverse_path, compressed=True)

    with open(output_json_read_count_file, 'w') as output_json_file:
        json.dump({'valid_reads_forward': valid_reads_forward,
                   'valid_reads_reverse': valid_reads_reverse},
                  output_json_file)

    if write_bam_files:
        bamfile_input = pysam.AlignmentFile(bamfile_input_path, "rb")
        bamfile_output_forward = pysam.AlignmentFile(output_bam_file_forward_path, "wb",
                                                     template=bamfile_input)
        bamfile_output_reverse = pysam.AlignmentFile(output_bam_file_reverse_path, "wb",
                                                     template=bamfile_input)
        for i, read in enumerate(bamfile_input):
            if i % 100_000 == 0:
                print(f"Writing reads: {i}", flush=True)
            if read.query_name in invalid_ids:
                continue
            if (read.is_read1 and read.is_forward) or (read.is_read2 and read.is_reverse):
                bamfile_output_forward.write(read)
            elif (read.is_read1 and read.is_reverse) or (read.is_read2 and read.is_forward):
                bamfile_output_reverse.write(read)
            else:
                assert False

        bamfile_output_forward.close()
        bamfile_output_reverse.close()
        bamfile_input.close()

        print("Sorting output .bam files")
        for bam_file_path in (output_bam_file_forward_path, output_bam_file_reverse_path):
            tmp_file_path = bam_file_path.parent / 'tmp_sorted.bam'
            pysam.sort("-o", str(tmp_file_path), str(bam_file_path))
            os.rename(tmp_file_path, bam_file_path)
            pysam.index(str(bam_file_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_file_input',
                        default='/cellfile/datapublic/jkoubele/compare_coverage/K002000093_54873Aligned.sortedByCoord.out.bam',
                        help='Folder containing .bam file to be processed.')
    # parser.add_argument('--bam_file_input',
    #                     default='/cellfile/datapublic/jkoubele/FLI_total_RNA/BAM_dedup_exact_position/no003-1_OA3/deduplicated.bam',
    #                     help='Folder containing .bam file to be processed.')
    # parser.add_argument('--bam_file_input',
    #                     default='/cellfile/datapublic/jkoubele/FLI_total_RNA/cov_example/chr1.bam',
    #                     help='Folder containing .bam file to be processed.')
    parser.add_argument('--output_folder',
                        default=None,
                        help='Output files will be written here.')
    args = parser.parse_args()
    extract_and_save_unique_pairs(Path(args.bam_file_input))
