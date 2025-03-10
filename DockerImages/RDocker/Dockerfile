# Base image
FROM r-base:latest

RUN apt-get update && apt-get install -y git build-essential

RUN git clone https://github.com/EPPIcenter/THEREALMcCOIL.git

WORKDIR /THEREALMcCOIL/categorical_method/

RUN rm *.so 

RUN rm *.o

RUN R CMD SHLIB McCOIL_categorical_code.c llfunction_het.c

RUN R -e "install.packages('readxl', repos='http://cran.r-project.org')"
