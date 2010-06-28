CC      =g++ -DGNU
CFLAGS  = -Wall -Wextra
DEBUG   = -g -O0 -DDEBUG -ggdb
POPCC   = popcc
METEOIODIR = ../src

LIBS      = -rdynamic -lstdc++ -ldl -lmeteoio
LIBS_POPC = -L../lib -lstdc++ -lmeteoiopopc

INCLUDE=-I. -I$(METEOIODIR)
INCLUDE_POPC=-I. -I$(METEOIODIR)

#####################
#	RULES
#####################
.cc.o: $*.cc $*.h
	$(CC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE)

%.o: %.cc
	$(CC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE)

%_par.o : %.cc
	$(POPCC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE_POPC) -o $@

 %.stub.o : %.ph
	$(POPCC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE_POPC) -o $@

#####################
#	TARGETS
#####################
seq: meteoio_demo

demo: meteoio_demo

par: meteoio_demo.popc

testDate: testDate.o 
	g++ -o testDate testDate.o $(LIBS)

testMeteo: testMeteo.o 
	g++ -o testMeteo testMeteo.o $(LIBS)

testStation: StationData.o testStation.o
	g++ -o testStation StationData.o testStation.o

meteoio_demo: meteoio_demo.o
	$(CC) $(DEBUG) -o $@ meteoio_demo.o ${LIBS}

meteoio_demo.popc: meteoio_demo_par.o 
	$(POPCC) $(DEBUG) -o $@ meteoio_demo_par.o ${LIBS_POPC}

testThomas: testThomas.o ../lib/libmeteoio.a
	cd ../; make 	
	$(CC) $(DEBUG) -rdynamic -o $@ testThomas.o ${LIBS} 

testwriting: testwriting.o
	$(CC) $(DEBUG) -o $@ testwriting.o ${LIBS}

testCoords: testCoords.o
	$(CC) $(DEBUG) -o $@ testCoords.o ${LIBS}

sniotest: sniotest.o
	$(CC) $(DEBUG) -o $@ sniotest.o ${LIBS}

sensorscope: sensorscope.o ../lib/libmeteoio.a
	$(CC) $(DEBUG) -rdynamic -o $@ sensorscope.o ${LIBS} 

testConfigReader: main.o ConfigReader.o slfutils.o slfexceptions.o
	$(CC) -o $@ main.o ConfigReader.o slfutils.o slfexceptions.o

hnw_interpolate: hnw_interpolate.o
	$(CC) $(DEBUG) -o $@ hnw_interpolate.o ${LIBS}

clean:
	rm -rf *~ *.o main testThomas meteoio_demo testDate testConfigReader testMeteo testStation meteoio_demo.popc *.out testArray sensorscope testwriting testCoords

