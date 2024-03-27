from ftplib import FTP
from pathlib import Path

from dotenv import dotenv_values
from tqdm import tqdm

if __name__ == "__main__":
    config = dotenv_values()
    ftp = FTP("genome.leibniz-fli.de")
    login_response = ftp.login(user=config['FLI_FTP_USER'], passwd=config['FLI_FTP_PASSWORD'])

    dataset_name = '20240219_866_YC'

    local_data_folder = Path('/cellfile/datapublic/jkoubele/FLI_total_RNA')
    remote_file_paths = ftp.nlst(f'/data/{dataset_name}')
    (local_data_folder / dataset_name).mkdir(exist_ok=True, parents=True)
    for file_path in tqdm(remote_file_paths, desc=f'Downloading {dataset_name}'):
        with open(local_data_folder / dataset_name / file_path.split('/')[-1], 'wb') as file:
            ftp.retrbinary(f'RETR {file_path}', file.write)
