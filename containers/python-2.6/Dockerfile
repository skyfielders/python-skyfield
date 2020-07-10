FROM mrupgrade/deadsnakes:2.6
RUN apt update
RUN apt install -y -y build-essential python2.6-dev
RUN pip install numpy==1.11.3
RUN pip install argparse mock pytz unittest2 pytest
RUN pip install skyfield
RUN echo 'assay --batch skyfield.tests' > /root/.bash_history
CMD cd skyfield/ci && /bin/bash
