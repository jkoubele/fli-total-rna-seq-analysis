import pandas as pd
from pathlib import Path
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import gseapy as gp

# %%
# Utility functions

L2FC_OF_50_PCT_INCREASE = np.log2(1.5)


def get_de_genes(df: pd.DataFrame,
                 direction: str = 'up',
                 pval_treshold=0.1,
                 lfold_min_size=L2FC_OF_50_PCT_INCREASE) -> set[str]:
    if direction == 'up':
        return set(df['gene_id'][(df['log2FoldChange'] > 0) &
                                 (np.abs(df['log2FoldChange']) > lfold_min_size) &
                                 (df['padj'] <= pval_treshold)])
    elif direction == 'down':
        return set(df['gene_id'][(df['log2FoldChange'] < 0) &
                                 (np.abs(df['log2FoldChange']) > lfold_min_size) &
                                 (df['padj'] <= pval_treshold)])
    else:
        assert False, "Direction must be either 'up' or 'down'"


# %%
# Load data  

gene_names = pd.read_csv('/cellfile/datapublic/jkoubele/STAR_2.7.11a/reference_genomes/GRCm38/gene_names.csv')
gene_name_by_id = {gene_id: gene_name for gene_id, gene_name in zip(gene_names['gene_id'], gene_names['gene_name'])}

folder_path = Path('/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/gene_DE_tables')

contrasts: dict[tuple[str]: pd.DataFrame] = {}
for file_path in folder_path.iterdir():
    compared_groups = file_path.stem.split('_vs_')
    contrasts[tuple(compared_groups)] = pd.read_csv(file_path).rename(columns={'Unnamed: 0': 'gene_id'})
    contrasts[tuple(compared_groups)]['gene_name'] = [gene_name_by_id[gene_id] for gene_id in
                                                      contrasts[tuple(compared_groups)]['gene_id']]

# %%
# Venn diagrams of DE genes (compared to AL group)        
dr_up = get_de_genes(contrasts[('DR', 'AL')])
nad_up = get_de_genes(contrasts[('NAD', 'AL')])
drnad_up = get_de_genes(contrasts[('DR_NAD', 'AL')])

dr_down = get_de_genes(contrasts[('DR', 'AL')], direction='down')
nad_down = get_de_genes(contrasts[('NAD', 'AL')], direction='down')
drnad_down = get_de_genes(contrasts[('DR_NAD', 'AL')], direction='down')

venn_diagram_folder = Path('/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/gene_venn_diagrams')

figure(figsize=(6, 5), dpi=180)

venn3([dr_up, nad_up, drnad_up], ('DR', 'NAD', 'DR_NAD'))
plt.title('Upregulated genes vs AL')
plt.savefig(venn_diagram_folder / 'upregulated.png')
plt.show()

figure(figsize=(6, 5), dpi=180)
venn3([dr_down, nad_down, drnad_down], ('DR', 'NAD', 'DR_NAD'))
plt.title('Downregulated genes vs AL')
plt.savefig(venn_diagram_folder / 'downregulated.png')
plt.show()

# %%
# Gene enrichment - load data

database_names = gp.get_library_name(organism='Mouse')

go_bp = gp.get_library(name='GO_Biological_Process_2023', organism='Mouse')
go_mf = gp.get_library(name='GO_Molecular_Function_2018', organism='Mouse')
go_cc = gp.get_library(name='GO_Cellular_Component_2023', organism='Mouse')
# %%
# Volcano plots + GSEA
volcano_plots_folder = Path('/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/volcano_plots')
enrichment_folder = Path('/cellfile/datapublic/jkoubele/FLI_total_RNA/DE_analysis/gene_enrichment')

custom_annotation_tresholds = {
    ('NAD', 'AL'): (4, 6),
    ('DR_NAD', 'DR'): (1, 6),
    ('DR', 'NAD'): (11, 7),
    ('DR_NAD', 'NAD'): (8.9, 5),
    ('DR_NAD', 'AL'): (17, 7),
    ('DR', 'AL'): (19, 6.6)}

for contrast_groups, df_contrast in contrasts.items():
    # contrast_groups = ('DR', 'AL')
    # df_contrast = contrasts[contrast_groups]

    mask_up = (df_contrast['padj'] < 0.1) & (df_contrast['log2FoldChange'] > 0) & (
            np.abs(df_contrast['log2FoldChange']) > L2FC_OF_50_PCT_INCREASE)
    mask_down = (df_contrast['padj'] < 0.1) & (df_contrast['log2FoldChange'] < 0) & (
            np.abs(df_contrast['log2FoldChange']) > L2FC_OF_50_PCT_INCREASE)
    mask_not_changed = ~ (mask_up | mask_down)

    annotation_mask = (-np.log10(df_contrast['padj']) > custom_annotation_tresholds[contrast_groups][0]) | (
            np.abs(df_contrast['log2FoldChange']) > custom_annotation_tresholds[contrast_groups][1])

    df_up = df_contrast[mask_up]
    df_down = df_contrast[mask_down]
    df_not_changed = df_contrast[mask_not_changed]

    figure(figsize=(6, 6), dpi=180)

    plt.scatter(x=df_up['log2FoldChange'],
                y=-np.log10(df_up['padj']),
                color='indianred',
                label='UP')

    plt.scatter(x=df_down['log2FoldChange'],
                y=-np.log10(df_down['padj']),
                color='deepskyblue',
                label='DOWN')

    plt.scatter(x=df_not_changed['log2FoldChange'],
                y=-np.log10(df_not_changed['padj']),
                color='lightgrey',
                label='NO CHANGE')
    plt.ylabel('-$\log_{10}$(p-value adj.)')
    plt.xlabel('L2FC')
    plt.title(f"{contrast_groups[0]} vs. {contrast_groups[1]}")

    plt.axhline(y=1, color='grey', linestyle='--', linewidth=1)
    plt.axvline(x=L2FC_OF_50_PCT_INCREASE, color='grey', linestyle='--', linewidth=1)
    plt.axvline(x=-L2FC_OF_50_PCT_INCREASE, color='grey', linestyle='--', linewidth=1)

    plt.legend()

    for x, y, gene_name in zip(df_contrast['log2FoldChange'][annotation_mask],
                               -np.log10(df_contrast['padj'][annotation_mask]),
                               df_contrast['gene_name'][annotation_mask]):
        plt.annotate(gene_name, (x, y))

    plt.savefig(volcano_plots_folder / f'{contrast_groups[0]}_vs_{contrast_groups[1]}.png')
    plt.show()

    background_genes = list(df_contrast['gene_name'].apply(lambda x: x.upper()))

    enrichment_gene_sets_mapping = {'gs_ind_0': 'GO BP',
                                    'gs_ind_1': 'GO MF',
                                    'gs_ind_2': 'GO CC'}

    upregulated_genes = list(df_up['gene_name'].apply(lambda x: x.upper()))
    if upregulated_genes:
        enrichment_upregulated = gp.enrich(gene_list=upregulated_genes,
                                           gene_sets=[go_bp, go_mf, go_cc],
                                           background=background_genes,
                                           outdir=None,
                                           verbose=True)
        enrichment_upregulated_df = enrichment_upregulated.results.sort_values(by='Adjusted P-value')
        enrichment_upregulated_df['Gene_set'] = [enrichment_gene_sets_mapping[x] for x in
                                                 enrichment_upregulated_df['Gene_set']]

        enrichment_upregulated_df.to_csv(
            enrichment_folder / f'{contrast_groups[0]}_vs_{contrast_groups[1]}_upregulated.csv')

    downregulated_genes = list(df_down['gene_name'].apply(lambda x: x.upper()))
    if downregulated_genes:
        enrichment_downregulated = gp.enrich(gene_list=downregulated_genes,
                                             gene_sets=[go_bp, go_mf, go_cc],
                                             background=background_genes,
                                             outdir=None,
                                             verbose=True)
        enrichment_downregulated_df = enrichment_downregulated.results.sort_values(by='Adjusted P-value')
        enrichment_downregulated_df['Gene_set'] = [enrichment_gene_sets_mapping[x] for x in
                                                   enrichment_downregulated_df['Gene_set']]

        enrichment_downregulated_df.to_csv(
            enrichment_folder / f'{contrast_groups[0]}_vs_{contrast_groups[1]}_downregulated.csv')
