FROM ubuntu:22.04

RUN apt-get update && apt-get install libboost-all-dev -y
RUN apt install libomp-dev -y
RUN apt install software-properties-common -y && add-apt-repository ppa:deadsnakes/ppa && apt update && apt install python3-pip -y
RUN apt-get install build-essential libssl-dev -y
RUN apt-get install cmake -y
RUN  apt install libnlopt-dev -y && apt install libnlopt-cxx-dev -y
WORKDIR /code
COPY ./ ./

RUN pip install CONET/python/conet-py

RUN mkdir build && cd build && cmake ../scicone/ && make
RUN pip install pyscicone/
ENTRYPOINT sleep infinity