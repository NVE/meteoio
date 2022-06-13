###
# Usage:
#--------------------
# Build this Dockerfile with MeteoIO.deb (Debian Release Build) being in the same folder:
# `$ sudo docker build -t meteoio_timeseries_web:latest .`
#
# Run the built image and mount a local folder for the MeteoIO job files:
# `$ sudo docker run --network=host --name=meteoio_timeseries_web -v /tmp/jobs:/jobs meteoio_timeseries_web:latest`
#
# Use `$ sudo docker ps` to list running Docker containers and `$ sudo docker stop <container>` to stop the conatiner.
# If you want to run a new conatiner, remove the old one with `$ sudo docker rm <container>`
#
###

# TODO: Use a more minimal image to build MeteoIO (e.g. alpine) and use that same image here
FROM gcc:latest

RUN apt-get update && apt-get -y install libzip-dev zipcmp zipmerge ziptool

WORKDIR /install

COPY ./MeteoIO.deb /install

RUN dpkg -i MeteoIO.deb && rm MeteoIO.deb

WORKDIR /usr/bin

# set up permissions for user `1001`
RUN chmod u+rx /usr/bin/meteoio_timeseries_web

EXPOSE 8080 8080

USER 1001

CMD ["./meteoio_timeseries_web", "-t 60", "-d /jobs"]
