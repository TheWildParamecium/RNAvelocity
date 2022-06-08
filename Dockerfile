FROM r-base:4.1.1
RUN mkdir -p /home/ubuntu/workspace && mkdir -p /home/ubuntu/workspace/lib/R/site-library && chmod -R 777 /home/ubuntu/workspace
WORKDIR /home/ubuntu/workspace
RUN apt-get update && apt-get install -y apt-utils && \
    apt-get install -y libgit2-dev libxml2-dev libfontconfig1-dev \
    libharfbuzz-dev libfribidi-dev libhdf5-dev libgsl-dev libfreetype6-dev \
    libpng-dev libtiff5-dev libjpeg-dev libgeos-dev libgdal-dev libcairo2-dev libxt-dev \
    libssl-dev libcurl4-openssl-dev curl tree jq htop texinfo vim man-db less
COPY install_packages.R /home/ubuntu/workspace/install_packages.R
RUN Rscript /home/ubuntu/workspace/install_packages.R
RUN  apt-get install -y fftw3 fftw3-dev pkg-config && \  
    R -e "install.packages('qqconf',dependencies=TRUE,repos='http://cran.rstudio.com/',lib='/home/ubuntu/workspace/lib/R/site-library')"  && \
    R -e "install.packages('Seurat',dependencies=TRUE,repos='http://cran.rstudio.com/',lib='/home/ubuntu/workspace/lib/R/site-library')"
COPY install_packages_v2.R /home/ubuntu/workspace/install_packages_v2.R
RUN Rscript /home/ubuntu/workspace/install_packages_v2.R && \
    rm /home/ubuntu/workspace/install_packages*
COPY ./Simulations/code/velosim_pseudo.r /home/ubuntu/workspace/velosim_pseudo.r
COPY ./Simulations/code/dyngen_pseudo.r /home/ubuntu/workspace/dyngen_pseudo.r
CMD ["/bin/bash"]