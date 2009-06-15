CC=g++ -DGNU
A3DROOT=../..
CFLAGS=-Wall -Wextra
PAROCC   = parocc

METEOIODIR	= ../src
FILTERDIR	= $(METEOIODIR)/filter

vpath    %.h     $(METEOIODIR) $(FILTERDIR) ./
vpath    %.c     $(METEOIODIR) $(FILTERDIR) ./
vpath    %.cc    $(METEOIODIR) $(FILTERDIR) ./

LIBS      =-rdynamic -lstdc++ -ldl -L../lib -lmeteoio -lfilter
LIBS_PAROC=-L../ -L../lib -L$(A3DROOT)/snowpack -L$(A3DROOT)/common -lstdc++ -lcommonparoc -lfilterparoc -lmeteoioparoc -lsnowpackparoc

DEBUG = -g -O0 -DDEBUG -ggdb
INCLUDE=-I. -I$(METEOIODIR) -I$(FILTERDIR)

INCLUDE_PAROC=-I. -I$(METEOIODIR) -I$(FILTERDIR) -I$(A3DROOT)/snowpack/ -I$(A3DROOT)/snowpack/snowpack_core/ -I$(A3DROOT)/common/ $(shell pkg-config --cflags libxml++-2.6) -I$(A3DROOT)/snowdrift

all: meteoio_demo

.cc.o: $*.cc $*.h
	$(CC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE)

%.o: %.cc
	$(CC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE)

testDate: Date.o testDate.o slfutils.o
	g++ -o testDate testDate.o Date.o slfutils.o

testMeteo: MeteoData.o Date.o testMeteo.o slfutils.o makeLib
	g++ -o testMeteo testMeteo.o MeteoData.o Date.o slfutils.o $(LIBS)

testStation: StationData.o testStation.o
	g++ -o testStation StationData.o testStation.o

#testFile: testFile.o makeLib
meteoio_demo: meteoio_demo.o
	$(CC) $(DEBUG) -o $@ meteoio_demo.o ${LIBS}

testThomas: testThomas.o 
	$(CC) $(DEBUG) -rdynamic -o $@ testThomas.o ${LIBS} 

testConfigReader: main.o ConfigReader.o slfutils.o slfexceptions.o
	$(CC) -o $@ main.o ConfigReader.o slfutils.o slfexceptions.o

all_paroc: meteoio_demo.paroc

%_par.o : %.cc
	$(PAROCC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE_PAROC) -o $@

 %.stub.o : %.ph
	$(PAROCC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE_PAROC) -o $@

meteoio_demo.paroc: meteoio_demo_par.o 
	$(PAROCC) $(DEBUG) -o $@ meteoio_demo_par.o ${LIBS_PAROC}

clean:
	rm -rf *~ *.o libslfio.a main testThomas meteoio_demo testDate testConfigReader testMeteo testStation meteoio_demo.paroc *.out testArray

