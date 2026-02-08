FROM ubuntu:20.04
COPY sources.list /etc/apt/sources.list
RUN apt update
RUN apt install -y -y curl python3-dev make
RUN curl -LsSf https://astral.sh/uv/0.10.0/install.sh \
    | UV_INSTALL_DIR="/usr/local/bin" sh

COPY requirements.txt /root/requirements.txt
RUN uv pip install --system -r /root/requirements.txt
RUN uv pip install --system jplephem sgp4
RUN ln -s /usr/bin/python3 /usr/local/bin/python

WORKDIR /src
ENV PYTHONDONTWRITEBYTECODE PLEASE
ENV PYTHONPATH /src
CMD ["/bin/bash"]
