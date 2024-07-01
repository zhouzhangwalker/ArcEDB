FROM ubuntu:20.04
#install build dependencies
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    g++-10 \
    apt-utils \
    ca-certificates \
    git \
    cmake \
    libgmp-dev \
    libomp-dev \
    libntl-dev

# Build ArcEDB
COPY . /ArcEDB
RUN mkdir /ArcEDB/build
WORKDIR /ArcEDB/build
RUN cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-10 -DSEAL_THROW_ON_TRANSPARENT_CIPHERTEXT=OFF -DSEAL_USE_INTEL_HEXL=ON && make -j4