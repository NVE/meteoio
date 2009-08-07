CC=g++ -DGNU
A3DROOT=../..
CFLAGS=-Wall -Wextra
POPCC   = popcc

METEOIODIR	= ../src
FILTERDIR	= $(METEOIODIR)/filter

vpath    %.h     $(METEOIODIR) $(FILTERDIR) ./
vpath    %.c     $(METEOIODIR) $(FILTERDIR) ./
vpath    %.cc    $(METEOIODIR) $(FILTERDIR) ./

LIBS      =-rdynamic -lstdc++ -ldl -L../lib -lmeteoio -lfilter
LIBS_POPC=-L../ -L../lib -L$(A3DROOT)/snowpack -L$(A3DROOT)/common -lstdc++ -lcommonpopc -lfilterpopc -lmeteoiopopc -lsnowpackpopc

DEBUG = -g -O0 -DDEBUG -ggdb
INCLUDE=-I. -I$(METEOIODIR) -I$(FILTERDIR)

INCLUDE_POPC=-I. -I$(METEOIODIR) -I$(FILTERDIR) -I$(A3DROOT)/snowpack/ -I$(A3DROOT)/snowpack/snowpack_core/ -I$(A3DROOT)/common/ $(shell pkg-config --cflags libxml++-2.6) -I$(A3DROOT)/snowdrift

all: testThomas

seq: meteoio_demo

.cc.o: $*.cc $*.h
	$(CC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE)

%.o: %.cc
	$(CC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE)

testDate: testDate.o ../lib/libmeteoio.a
	cd ../; make
	g++ -o testDate testDate.o $(LIBS)

testMeteo: MeteoData.o Date.o testMeteo.o slfutils.o makeLib
	g++ -o testMeteo testMeteo.o MeteoData.o Date.o slfutils.o $(LIBS)

testStation: StationData.o testStation.o
	g++ -o testStation StationData.o testStation.o

#testFile: testFile.o makeLib
meteoio_demo: meteoio_demo.o
	$(CC) $(DEBUG) -o $@ meteoio_demo.o ${LIBS}

testThomas: testThomas.o ../lib/libmeteoio.a
	cd ../; make 	
	$(CC) $(DEBUG) -rdynamic -o $@ testThomas.o ${LIBS} 

sensorscope: sensorscope.o ../lib/libmeteoio.a
	$(CC) $(DEBUG) -rdynamic -o $@ sensorscope.o ${LIBS} 

testConfigReader: main.o ConfigReader.o slfutils.o slfexceptions.o
	$(CC) -o $@ main.o ConfigReader.o slfutils.o slfexceptions.o

par: meteoio_demo.popc

%_par.o : %.cc
	$(POPCC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE_POPC) -o $@

 %.stub.o : %.ph
	$(POPCC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE_POPC) -o $@

meteoio_demo.popc: meteoio_demo_par.o 
	$(POPCC) $(DEBUG) -o $@ meteoio_demo_par.o ${LIBS_POPC}

clean:
	rm -rf *~ *.o libslfio.a main testThomas meteoio_demo testDate testConfigReader testMeteo testStation meteoio_demo.popc *.out testArray sensorscope

