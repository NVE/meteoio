CC      = g++ -DGNU
CFLAGS  = -Wall -Wextra
DEBUG   = -g -O0 -DDEBUG -ggdb
METEOIODIR = ../../src

LIBS      = -rdynamic -lstdc++ -ldl -L../../lib -lmeteoio

INCLUDE=-I. -I$(METEOIODIR)
INCLUDE_POPC=-I. -I$(METEOIODIR)

#####################
#	RULES
#####################
.cc.o: $*.cc $*.h
	$(CC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE)

%.o: %.cc
	$(CC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE)

#####################
#	TARGETS
#####################
all: meteo_reading dem_reading

meteo_reading: meteo_reading.o
	$(CC) $(DEBUG) -o $@ meteo_reading.o ${LIBS}

dem_reading: dem_reading.o
	$(CC) $(DEBUG) -o $@ dem_reading.o ${LIBS}


clean:
	rm -rf *~ *.o  *.out meteo_reading dem_reading

