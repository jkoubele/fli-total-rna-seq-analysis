# fli-total-rna-seq-analysis

## Setup

- To install required Python packages, please run ```pip install -r requirements.txt```.
- If you wish to download the data from the FLI FTP server, please copy the file ```.env_example``` and
  store it as a new file ```.env``` (which is not meant to be commited to git). Then, fill the required credentials in
  the ```.env``` file. (The download script
  ```download_data.py``` is using these credentials.)