FROM ubuntu:bionic

# Set the working directory
WORKDIR /reproducibility

# Copy the directory contents into the container 
COPY . /reproducibility

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get -y install cmake g++ python3-pip libbz2-dev libboost-all-dev wget bash \
  && pip3 install --trusted-host pypi.python.org -r requirements.txt

# Compile the software
RUN mkdir build \
  && cd build \
  && cmake .. \
  && make

# Download and compile MCL
RUN mkdir mcl-build \
  && cd mcl-build \
  && wget https://micans.org/mcl/src/mcl-14-137.tar.gz \
  && tar xvf mcl-14-137.tar.gz \
  && cd mcl-14-137 \
  && ./configure --prefix=/reproducibility/mcl-build \
  && make install

# Make the directory to hold the results
RUN mkdir /results

ENV PATH="/reproducibility/mcl-build/bin:${PATH}"

CMD cd /reproducibility/Reproducibility && bash run.sh && cp -r * /results

