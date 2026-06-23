# Base image with R 4.5.1
FROM rocker/shiny:4.5.1

# General updates and system dependencies
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y build-essential gfortran git wget curl \
    libxml2-dev libmagick++-dev libssl-dev libharfbuzz-dev libfribidi-dev \
    libcurl4-openssl-dev libgsl-dev libgit2-dev libssh2-1-dev libudunits2-dev \
    libpng-dev libjpeg-dev libfreetype6-dev libtiff-dev liblapack-dev libblas-dev \
    pandoc libopenblas-dev libreadline-dev libncurses-dev zlib1g-dev \
    libhdf5-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install R packages directly
RUN Rscript -e 'install.packages(c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "magrittr", "ggplot2", "ggrepel", "hdf5r", "ggdendro", "gridExtra"), dependencies = TRUE)'

# Copy the app files (scripts, data, etc.)
COPY shinyAppMulti/ /srv/shiny-server/

# Ensure that the expected user is present in the container
RUN if id shiny &>/dev/null && [ "$(id -u shiny)" -ne 999 ]; then \
        userdel -r shiny; \
        id -u 999 &>/dev/null && userdel -r "$(id -un 999)"; \
    fi; \
    useradd -u 999 -m -s /bin/bash shiny; \
    chown -R shiny:shiny /srv/shiny-server/ /var/lib/shiny-server/ /var/log/shiny-server/

# Other settings
USER shiny
EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
