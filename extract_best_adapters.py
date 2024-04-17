import pandas as pd
from pathlib import Path
import argparse


def aggregate_adapters(folder_path: Path, output_file_name='aggregated_adapters.txt') -> None:
    adapters_1: set[str] = set()
    adapters_2: set[str] = set()
    for file in folder_path.iterdir():
        if output_file_name in file.name:
            continue  # we write the output to the same directory, so we skip such file if it already exists
        data = pd.read_csv(file, delimiter='\t').squeeze()
        adapters_1.add(data['best_adapter1'])
        adapters_2.add(data['best_adapter2'])

    error_message = "DETECTED ADAPTERS ARE INCONSISTENT BETWEEN SAMPLES!"
    best_adapter_1 = adapters_1.pop() if len(adapters_1) == 1 else error_message
    best_adapter_2 = adapters_2.pop() if len(adapters_2) == 1 else error_message
    with open(folder_path / output_file_name, 'w') as out_file:
        out_file.write(f"{best_adapter_1}\n{best_adapter_2}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--detected_adapters_folder',
                        default='/cellfile/datapublic/jkoubele/FLI_total_RNA/detected_adapters',
                        help='Folder containing output Atria adapters detection.')
    args = parser.parse_args()
    aggregate_adapters(folder_path=Path(args.detected_adapters_folder))
