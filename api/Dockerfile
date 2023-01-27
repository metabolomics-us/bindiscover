FROM continuumio/miniconda3

COPY binvestigate_api.yml .

RUN conda env create -f binvestigate_api.yml

RUN mkdir -p /newer_api/assets

RUN mkdir -p /newer_api/assets

RUN mkdir -p /newer_api/pages

RUN mkdir -p /newer_datasets/

COPY ./newer_api/*.py /newer_api/

COPY ./newer_datasets/* /newer_datasets/

WORKDIR /newer_api

EXPOSE 4999

SHELL ["conda", "run", "-n", "binvestigate_api", "/bin/bash", "-c"]

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "binvestigate_api", "python", "./parent_api.py"]


