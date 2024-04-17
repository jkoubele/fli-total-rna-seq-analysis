FROM debian:latest
RUN apt update && apt upgrade -y
RUN apt install curl -y
RUN curl -L -O https://github.com/cihga39871/Atria/releases/download/v4.0.3/atria-4.0.3-linux-ubuntu22.tar.gz
RUN tar -zxf atria-4.0.3-linux-ubuntu22.tar.gz
RUN apt install pigz pbzip2 -y