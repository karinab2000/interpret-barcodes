# Base image
FROM python:3.11-slim

# Set environment variables
ENV PYTHONUNBUFFERED=1

# Set the working directory
WORKDIR /app

# Add in files
COPY barcode_pos.txt regression_models.pkl /usr/local/bin/

# Copy your application files into the container
COPY . /app

# Install dependencies
RUN pip install --no-cache-dir \
    pandas numpy openpyxl xlrd>=2.0.1 xlsxwriter \
    seaborn matplotlib scikit-learn scipy statsmodels networkx patsy \
    dill numpyencoder

RUN apt-get update && apt-get install -y \
    gcc g++ \
    && rm -rf /var/lib/apt/lists/*