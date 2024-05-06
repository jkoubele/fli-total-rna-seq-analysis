from enum import StrEnum
from pathlib import Path
from typing import Optional

import pandas as pd
from tqdm import tqdm

project_folder = Path('/cellfile/datapublic/jkoubele/FLI_total_RNA')
gene_counts_folder = project_folder / 'gene_counts'
gene_counts_htseq_folder = project_folder / 'gene_counts_htseq'
aggregate_counts_folder = project_folder / 'aggregate_counts'


def aggregate_feature_counts() -> None:
    """
    Aggregate all outputs from featureCounts into a single file.
    """
    sample_to_counts: dict[str, pd.Series] = {}
    template_df: Optional[pd.DataFrame] = None
    for sample_folder in tqdm(gene_counts_folder.iterdir(), desc='Aggregating feature counts to a single file'):
        sample_df = pd.read_csv(sample_folder / 'feature_counts.tsv', delimiter='\t', header=1)
        sample_to_counts[sample_folder.name] = sample_df['/bam_folder/deduplicated.bam']
        if template_df is None:
            template_df = sample_df.drop(columns='/bam_folder/deduplicated.bam')
        else:
            assert template_df.equals(sample_df.drop(columns='/bam_folder/deduplicated.bam'))
    for sample_name, counts in sample_to_counts.items():
        template_df[sample_name] = counts
    aggregate_counts_folder.mkdir(exist_ok=True)
    template_df.to_csv(aggregate_counts_folder / 'aggregate_feature_counts.tsv', sep='\t')


def aggregate_htseq_counts() -> None:
    """
    Aggregate all outputs from htseq-count into a single file.
    """
    sample_to_counts: dict[str, pd.Series] = {}
    template_df: Optional[pd.DataFrame] = None
    for sample_file in tqdm(gene_counts_htseq_folder.iterdir(),
                            desc='Aggregating htseq-count outputs to a single file'):
        sample_df = pd.read_csv(sample_file, delimiter='\t', names=['gene_id', 'count'])
        sample_to_counts[sample_file.stem] = sample_df['count']
        if template_df is None:
            template_df = sample_df.drop(columns='count')
        else:
            assert template_df.equals(sample_df.drop(columns='count'))
    for sample_name, counts in sample_to_counts.items():
        template_df[sample_name] = counts
    aggregate_counts_folder.mkdir(exist_ok=True)
    template_df.to_csv(aggregate_counts_folder / 'aggregate_counts_htseq.tsv', sep='\t')


class InterventionGroups(StrEnum):
    AL = 'AL'
    DR = 'DR'
    NAD = 'NAD'
    DR_NAD = 'DR_NAD'


def sample_name_to_intervention_group(sample_name: str) -> str:
    if 'OA' in sample_name:
        return InterventionGroups.AL.value
    elif 'ND' in sample_name:
        return InterventionGroups.DR_NAD.value
    elif 'N' in sample_name:
        return InterventionGroups.NAD.value
    return InterventionGroups.DR.value


def prepare_annontation_file() -> None:
    df_feature_counts = pd.read_csv(aggregate_counts_folder / 'aggregate_feature_counts.tsv', sep='\t')
    df_feature_counts = df_feature_counts[[column for column in df_feature_counts.columns if column.startswith('no')]]
    annotation_df = pd.DataFrame(data={'sample_name': df_feature_counts.columns,
                                       'group': [sample_name_to_intervention_group(sample) for sample in
                                                 df_feature_counts.columns]})
    annotation_df = annotation_df.set_index('sample_name')
    annotation_df.to_csv(aggregate_counts_folder / 'annotation.tsv', sep='\t')


if __name__ == "__main__":
    aggregate_feature_counts()
    aggregate_htseq_counts()
    prepare_annontation_file()
