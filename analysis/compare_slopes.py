import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from dataclasses_json import DataClassJsonMixin
from scipy import stats


@dataclass
class Group(DataClassJsonMixin):
    name: str
    sample_names: list[str]


@dataclass
class Experiment(DataClassJsonMixin):
    group_control: Group
    group_intervention: Group
    p_value_2_tail: float
    p_value_greater: float
    p_value_less: float


def load_valid_slopes(input_folder: Path, clip_by_zero=False) -> pd.DataFrame:
    input_df: dict[str, pd.DataFrame] = {}
    for file in input_folder.iterdir():
        if not file.suffix == '.tsv':
            continue
        input_df[file.stem] = pd.read_csv(file, sep='\t')
    df_slopes = pd.DataFrame(data={name: df['slope'] for name, df in input_df.items()})
    df_slopes = df_slopes[np.all(df_slopes.isna() == False, axis=1)]
    if clip_by_zero:
        df_slopes[df_slopes > 0] = 0
    df_slopes = df_slopes[np.all(df_slopes <= 0, axis=1)]
    return df_slopes


group_al = Group(name='AL',
                 sample_names=['no003-1_OA3', 'no019-0_OA4', 'no020-0_OA5', 'no021-0_OA6'])
group_dr = Group(name='DR',
                 sample_names=['no024-0_OD6', 'no005-0_OD2', 'no004-0_OD1', 'no022-0_OD4'])
group_nad = Group(name='NAD',
                  sample_names=['no025-0_ON4', 'no034-0_ON3', 'no026-0_ON5', 'no008-0_ON2'])
group_combined = Group(name='DR+NAD',
                       sample_names=['no030-0_OND6', 'no028-0_OND4', 'no010-0_OND1', 'no011-0_OND2'])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder', default="/cellfile/datapublic/jkoubele/FLI_total_RNA/intronic_slopes")
    parser.add_argument('--output_folder', default="/cellfile/datapublic/jkoubele/FLI_total_RNA/slopes_comparison")
    args = parser.parse_args()

    input_data_folder = Path(args.input_folder)
    output_folder = Path(args.output_folder)
    output_folder.mkdir(exist_ok=True, parents=True)

    df_slopes = load_valid_slopes(input_data_folder)
    for group_control, group_intervention in [(group_al, group_dr),
                                              (group_al, group_nad),
                                              (group_al, group_combined),
                                              (group_nad, group_combined),
                                              (group_dr, group_combined),
                                              (group_nad, group_dr)]:
        avg_slope_control = df_slopes[group_control.sample_names].mean(axis=1)
        avg_slope_intervention = df_slopes[group_intervention.sample_names].mean(axis=1)
        experiment = Experiment(group_control=group_control,
                                group_intervention=group_intervention,
                                p_value_2_tail=stats.wilcoxon(avg_slope_control, avg_slope_intervention,
                                                              alternative='two-sided').pvalue,
                                p_value_greater=stats.wilcoxon(avg_slope_control, avg_slope_intervention,
                                                               alternative='greater').pvalue,
                                p_value_less=stats.wilcoxon(avg_slope_control, avg_slope_intervention,
                                                            alternative='less').pvalue)
        with open(output_folder / f"{experiment.group_control.name}_vs_{experiment.group_intervention.name}.json",
                  'w') as output_file:
            output_file.write(experiment.to_json())

        print(experiment)
        print(50 * '-')
