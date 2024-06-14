FROM python:3.11
COPY . .
RUN apt update && apt upgrade -y

# Install bedtools
RUN curl -L -O https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz
RUN tar -xvzf bedtools-2.31.1.tar.gz
RUN cd /bedtools2 && make
RUN echo 'PATH=$PATH:/bedtools2/bin/' >> ~/.bashrc

RUN pip install -r requirements.txt

