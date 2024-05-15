from enum import StrEnum
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from sklearn.decomposition import PCA


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


colors_by_intervention_group = {InterventionGroups.AL.value: 'red',
                                InterventionGroups.DR.value: 'blue',
                                InterventionGroups.NAD.value: 'yellow',
                                InterventionGroups.DR_NAD.value: 'green'}

project_folder = Path('/cellfile/datapublic/jkoubele/FLI_total_RNA')
aggregate_counts_folder = project_folder / 'aggregate_counts'

df_feature_counts = pd.read_csv(aggregate_counts_folder / 'aggregate_feature_counts.tsv', sep='\t')
df_htseq_counts = pd.read_csv(aggregate_counts_folder / 'aggregate_counts_htseq.tsv', sep='\t')

df_feature_counts = df_feature_counts.rename(columns={'Geneid': 'gene_id'})
df_feature_counts = df_htseq_counts

df_feature_counts = df_feature_counts.set_index('gene_id')
df_feature_counts = df_feature_counts[[column for column in df_feature_counts.columns if column.startswith('no')]]

keep_mask = [gene_id.startswith('ENSMUSG') for gene_id in df_feature_counts.index]
df_feature_counts = df_feature_counts[keep_mask]

# drop rows with all zeros
df_feature_counts = df_feature_counts.loc[(df_feature_counts != 0).any(axis=1)]

intervention_groups = [sample_name_to_intervention_group(sample_name) for sample_name in df_feature_counts.columns]
intervention_group_number = [list(InterventionGroups).index(x) for x in intervention_groups]

X = df_feature_counts.copy().values

num_reads = X.sum(axis=0)

X_transformed = X / num_reads * 1e6
X_transformed = X_transformed.T
X_transformed = np.log1p(X_transformed)

X_transformed -= np.mean(X_transformed, axis=0)

std = np.std(X_transformed, axis=0)

# X_transformed /= std

pca = PCA().fit(X_transformed)
print(f"Explained variance \n PC1: {pca.explained_variance_ratio_[0]} \n PC2: {pca.explained_variance_ratio_[1]}")

X_pca = pca.transform(X_transformed)

colors = ListedColormap(['r', 'b', 'g', 'yellow'])
for group in set(intervention_groups):
    group_mask = [x == group for x in intervention_groups]
    scatter = plt.scatter(x=X_pca[:, 0][group_mask], y=X_pca[:, 1][group_mask], c=colors_by_intervention_group[group],
                          label=group)
# plt.legend(handles=scatter.legend_elements(), labels=intervention_groups)
plt.xlabel(f"PC 1 ({round(pca.explained_variance_ratio_[0] * 100, 1)}%)")
plt.ylabel(f"PC 2 ({round(pca.explained_variance_ratio_[1] * 100, 1)}%)")
plt.title('Log transform only')
plt.legend()
plt.show()
