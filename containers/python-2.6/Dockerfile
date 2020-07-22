FROM mrupgrade/deadsnakes:2.6
RUN apt update
RUN apt install -y -y build-essential python2.6-dev
RUN pip install numpy==1.11.3
RUN pip install argparse certifi jplephem mock pytz sgp4 unittest2
RUN pip install https://github.com/brandon-rhodes/assay/archive/master.zip
RUN echo 'PYTHONPATH=.. assay --batch skyfield.tests' > /root/.bash_history
CMD cd skyfield/ci && /bin/bash
