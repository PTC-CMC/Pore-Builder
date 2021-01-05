FROM continuumio/miniconda3:4.8.2-alpine AS builder

ENV PATH /opt/conda/bin:$PATH

USER root

RUN apk --update-cache upgrade && \
    apk add --no-cache bash vim nano && \
    apk add --no-cache --virtual .build_deps cmake make git git-lfs && \
	mkdir /software && cd /software && \
	conda update conda -y && \
	. /opt/conda/etc/profile.d/conda.sh && \
    conda create --name pore37 mbuild mdtraj parmed openbabel \
    -c conda-forge -c anaconda && \
	echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/anaconda/.bashrc && \
	echo "conda activate base" >> /home/anaconda/.bashrc && \
	echo "conda activate pore37" >> /home/anaconda/.bashrc && \
	conda activate pore37 && \
	cd /software && \
    git clone https://github.com/rmatsum836/Pore-Builder && \
    cd Pore-Builder && \
    python setup.py install && cd ../ && \
	conda clean -afy && apk del .build_deps && \
    rm -rf /var/cache/apk/*

WORKDIR /workspace

COPY entrypoint.sh /entrypoint.sh
RUN chmod a+x /entrypoint.sh && \
    chown -R anaconda /workspace && \
    chmod 755 /workspace

USER anaconda

ENTRYPOINT ["/entrypoint.sh"]
CMD ["none"]
