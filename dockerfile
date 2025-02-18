FROM python:3.9

WORKDIR /app
COPY analyze.py /app/
COPY requirements.txt /app/

RUN pip install -r requirements.txt

ENTRYPOINT ["/bin/bash"]

