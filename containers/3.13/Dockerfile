FROM python:3.13-slim
WORKDIR /work

RUN pip install https://github.com/brandon-rhodes/assay/archive/master.zip
RUN pip install certifi jplephem numpy pandas sgp4

CMD ["bash", "-c", "PYTHONPATH=$PWD ./test-code.sh"]
