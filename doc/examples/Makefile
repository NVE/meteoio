CC      = g++ -DGNU
CFLAGS  = -Wall -Wextra
DEBUG   = -g #-O3 # -DDEBUG -ggdb
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
all: meteo_reading dem_reading coordinates 2D_interpolations time

meteo_reading: meteo_reading.o
	$(CC) $(DEBUG) -o $@ meteo_reading.o ${LIBS}
dem_reading: dem_reading.o
	$(CC) $(DEBUG) -o $@ dem_reading.o ${LIBS}
coordinates: coordinates.o
	$(CC) $(DEBUG) -o $@ coordinates.o ${LIBS}
2D_interpolations: 2D_interpolations.o
	$(CC) $(DEBUG) -o $@ 2D_interpolations.o ${LIBS}
time: time.o
	$(CC) $(DEBUG) -o $@ time.o ${LIBS}

clean:
	rm -rf *~ *.o  *.out meteo_reading dem_reading coordinates 2D_interpolations time

