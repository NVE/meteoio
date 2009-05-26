CC=g++ -DGNU
A3DROOT=../..
CFLAGS=-Wall -Wextra
PAROCC   = parocc

SLFIODIR	= ../src
FILTERDIR	= $(SLFIODIR)/filter

vpath    %.h     $(SLFIODIR) $(FILTERDIR) ./
vpath    %.c     $(SLFIODIR) $(FILTERDIR) ./
vpath    %.cc    $(SLFIODIR) $(FILTERDIR) ./

LIBS      =-L../ -L../../snowpack -rdynamic -lstdc++ -lsnowpack -ldl -L../lib -lmeteoio -lfilter 
LIBS_PAROC=-L../ -L../filter -L../../snowpack -L../../common -lstdc++ -lcommonparoc -lfilterparoc -lmeteoioparoc -lsnowpackparoc

DEBUG = -g -O0 -DDEBUG -ggdb
INCLUDE=-I. -I$(SLFIODIR) -I$(FILTERDIR) -I$(A3DROOT)/snowpack/ -I$(A3DROOT)/snowpack/snowpack_core/ -I$(A3DROOT)/common/ 

INCLUDE_PAROC=-I. -I$(SLFIODIR) -I$(FILTERDIR) -I$(A3DROOT)/snowpack/ -I$(A3DROOT)/snowpack/snowpack_core/ -I$(A3DROOT)/common/ $(shell pkg-config --cflags libxml++-2.6) -I$(A3DROOT)/snowdrift

all: testFile

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
testFile: testFile.o
	$(CC) $(DEBUG) -o $@ testFile.o ${LIBS}

testThomas: testThomas.o 
	$(CC) $(DEBUG) -rdynamic -o $@ testThomas.o ${LIBS} 

testConfigReader: main.o ConfigReader.o slfutils.o slfexceptions.o
	$(CC) -o $@ main.o ConfigReader.o slfutils.o slfexceptions.o

all_paroc: testFile.paroc

%_par.o : %.cc
	$(PAROCC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE_PAROC) -o $@

 %.stub.o : %.ph
	$(PAROCC) $(DEBUG) $(CFLAGS) -c $< $(INCLUDE_PAROC) -o $@

testFile.paroc: testFile_par.o 
	$(PAROCC) $(DEBUG) -o $@ testFile_par.o ${LIBS_PAROC}

clean:
	rm -rf *~ *.o libslfio.a main testThomas testFile testDate testConfigReader testMeteo testStation testFile.paroc *.out testArray

