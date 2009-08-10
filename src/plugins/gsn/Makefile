CPP=g++
#INCLUDE=-I /usr/include/gsoap/
INCLUDE=
OBJECTS=test.o soapA3DWebServiceSoap12BindingProxy.o soapC.o
#OBJECTS=test.o
LIBS=-lgsoap++

.SUFFIXES: .o .cpp .cc

.cpp.o: $*.cpp $*.h
	$(CPP) -c $< $(INCLUDE) -o $@

.cc.o: $*.cc $*.h
	$(CPP) -c $< $(INCLUDE) -o $@

all: wsdl2cpp client

soapA3DWebServiceSoap12BindingProxy.o:
	$(CPP) -c soapA3DWebServiceSoap12BindingProxy.cpp $(INCLUDE) -o $@

soapC.o:
	$(CPP) -c soapC.cpp $(INCLUDE) -o $@

client: $(OBJECTS)
	$(CPP) -o $@ $(OBJECTS) $(LIBS)

wsdl2cpp:
	wsdl2h -g A3DWebService.wsdl
	soapcpp2 -2 -C -L -x -w -i A3DWebService.h $(INCLUDE)
	rm -f *Soap11*

clean:
	rm -f *~ *.o soap* *.nsmap client A3DWebService.h

