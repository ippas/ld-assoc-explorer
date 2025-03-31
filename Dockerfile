FROM python:3.6.9-slim-stretch

WORKDIR /usr/src/app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

#RUN useradd app
#USER app

ENTRYPOINT [ "./analysis/run-analysis.sh"]