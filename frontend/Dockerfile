FROM continuumio/miniconda3

COPY binvestigate_gui_2.yml .

RUN conda env create -f binvestigate_gui_2.yml \
    &&  mkdir -p /newer_frontend/assets \
    &&  mkdir -p /newer_frontend/pages \
    &&  mkdir -p /newer_datasets/ \
    &&  conda clean -afy


COPY ./newer_frontend/app.py /newer_frontend/

COPY ./newer_frontend/pages/*.py /newer_frontend/pages/

COPY ./newer_frontend/assets/* /newer_frontend/assets/

COPY ./newer_datasets/* /newer_datasets/

WORKDIR /newer_frontend

EXPOSE 8050

SHELL ["conda", "run", "-n", "binvestigate_gui", "/bin/bash", "-c"]

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "binvestigate_gui_2", "python", "./app.py"]


