#
# Makefile for the METEOIO Library
#
#################################################
# USER OPTIONS: MAKE CHANGES IN THIS SECTION ONLY
#################################################

####### COMPILERS AND OPTIONS
#Destination system: either zeus or grid
DEST	 = grid
#DEST	 = zeus

# build mode: release or debug
MODE	 = release
MODE	 = debug

# specific plugins to build: yes or no
BOSCHUNGIO	= no
IMISIO		= no

#SVNREV_GEN	= $(shell main/version.sh)
#SVNREV		= $(eval SVNREV := $(SVNREV_GEN))$(SVNREV)

ifeq ($(DEST),grid)		#for grid
	CXX      = g++
	#CXX	 = colorgcc
	LINKER   = g++ -DGNU
	PAROCC   = parocc
	#ARCH     = -march=prescott
	ARCH     = -march=pentium4
	OPTIM    = -fomit-frame-pointer -O3 #-fdata-sections
else				#for Zeus
        CXX      = g++
        #CXX     = pathcc -DGNU
        LINKER   = g++
        PAROCC   = parocc
        ARCH     = -march=x86-64
        OPTIM    = #-fomit-frame-pointer #-O3 -fdata-sections
        #OPTIM   = -O2 -Wstrict-aliasing #-fno-strict-aliasing
endif

DEBUG    = -g -O0 -D__DEBUG #-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC #-dr 
WARNING  = -Wall -Wextra #-Weffc++ #-pedantic-errors #-Werror

#DIRECTORIES
TARGET		= meteoIO
TARGET_PAROC	= meteoIO_paroc

SRCDIR		= ./src
LIBDIR		= ./lib
FILTERDIR	= $(SRCDIR)/filter

#################################################
# END USER OPTIONS
#################################################

ifeq ($(MODE),debug)		#for debug
	CCFLAGS  = -combine -pipe $(ARCH) $(DEBUG) $(WARNING)
else				#for release
	CCFLAGS  = -combine -pipe $(ARCH) $(OPTIM) $(WARNING)
endif
LIBS		= -lc
LDFLAGS_SEQ	=  -L$(LIBDIR) -lmeteoIO -lfilter
LD_PAROC	=  -L$(LIBDIR) -lmeteoIO -lfilter
LDFLAGS_PAROC	= $(LD_PAROC)
LDFLAGS		= $(LIBS)
INCLUDE		= -I$(SRCDIR) -I$(FILTERDIR)

INCLUDE_ORA = -I/software/oracle/client_1/rdbms/public/
LIBS_ORA = -L/software/oracle/client_1/lib -L/usr/lib/
LDFLAGS_ORA = -Xlinker -zmuldefs -locci -lclntsh -lstdc++ 
#CCFLAGS += -DSVNREV="\"$(SVNREV)\""


######## Sources, objects, headers
METEOIO_OBJ = 	$(SRCDIR)/MeteoData.o \
		$(SRCDIR)/StationData.o \
		$(SRCDIR)/IOHandler.o \
		$(SRCDIR)/IOInterface.o \
		$(SRCDIR)/DynamicLibrary.o \
		$(SRCDIR)/A3DIO.o \
		$(SRCDIR)/ConfigReader.o \
		$(SRCDIR)/Date_IO.o \
		$(SRCDIR)/Grid2DObject.o \
		$(SRCDIR)/IOExceptions.o \
		$(SRCDIR)/IOUtils.o \
		$(SRCDIR)/libinterpol2D.o \
		$(SRCDIR)/libinterpol1D.o \
		$(SRCDIR)/Meteo1DResampler.o \
		$(SRCDIR)/Meteo2DInterpolator.o \
		$(SRCDIR)/MeteoBuffer.o \
		$(SRCDIR)/Meteo1DResampler.o \
		$(SRCDIR)/LegacyIO.o 

METEOIO_OBJ_PAROC =  $(SRCDIR)/IOInterface_par.o \
		$(SRCDIR)/LegacyIO_par.o \
		$(SRCDIR)/IOInterface_par.o \
		$(SRCDIR)/IOHandler_par.o \
		$(SRCDIR)/IOInterface.stub.o \
		$(SRCDIR)/LegacyIO.stub.o \
		$(SRCDIR)/MeteoData_par.o \
		$(SRCDIR)/StationData_par.o \
		$(SRCDIR)/ConfigReader_par.o \
		$(SRCDIR)/Date_IO_par.o \
		$(SRCDIR)/Grid2DObject_par.o \
		$(SRCDIR)/IOExceptions_par.o \
		$(SRCDIR)/IOUtils_par.o \
		$(SRCDIR)/libinterpol2D_par.o \
		$(SRCDIR)/libinterpol1D_par.o \
		$(SRCDIR)/Meteo1DResampler_par.o \
		$(SRCDIR)/Meteo2DInterpolator_par.o \
		$(SRCDIR)/MeteoBuffer_par.o \
		$(SRCDIR)/Meteo1DResampler_par.o \
		$(SRCDIR)/DynamicLibrary_par.o 

FILTER_OBJ = 	$(FILTERDIR)/FilterBase.o \
		$(FILTERDIR)/FilterBase1Stn.o \
		$(FILTERDIR)/FilterValue.o \
		$(FILTERDIR)/MinValue.o \
		$(FILTERDIR)/MaxValue.o \
		$(FILTERDIR)/MinMaxValue.o \
		$(FILTERDIR)/MaxChangeRate.o \
		$(FILTERDIR)/NoObservedChange.o \
		$(FILTERDIR)/FilterFacade.o 

FILTER_OBJ_PAROC  = 	$(FILTERDIR)/FilterBase_par.o \
		$(FILTERDIR)/FilterBase1Stn_par.o \
		$(FILTERDIR)/FilterValue_par.o \
		$(FILTERDIR)/MinValue_par.o \
		$(FILTERDIR)/MaxValue_par.o \
		$(FILTERDIR)/MinMaxValue_par.o \
		$(FILTERDIR)/MaxChangeRate_par.o \
		$(FILTERDIR)/NoObservedChange_par.o \
		$(FILTERDIR)/FilterFacade_par.o 


####### Build rules SEQ

all: seq 
	@printf "*** MeteoIO compiled for \033[36m%s\033[0m\n" $(DEST)

help:
	@printf "MeteoIO Makefile targets for \033[36m%s\033[0m:\n" $(DEST)
	@printf " \033[36mseq\033[0m \n"
	@printf " \033[36mparoc\033[0m \n"
	@printf " \033[36minstall\033[0m \n"
	@printf " \033[36mclean\033[0m \n"
	@printf " \033[36mdistclean\033[0m \n"
	@printf " \033[36mdocumentation\033[0m\n"

seq: $(LIBDIR)/libmeteoIO.a build_dynamiclibs build_staticlibs

build_staticlibs: $(LIBDIR)/libfilter.a

build_dynamiclibs: $(LIBDIR)/libBoschungIO.so $(LIBDIR)/libImisIO.so

############## PAROC ##############
paroc: meteoIO_lib_paroc meteoIO_module_paroc filter_lib_paroc

filter_lib_paroc: libfilterparoc.a

meteoIO_lib_paroc: libmeteoIOparoc.a 

meteoIO_module_paroc: meteoIO.module
##############  END  ##############

clean:
	rm -f $(SRCDIR)/*~ $(LIBDIR)/*.a $(SRCDIR)/*.o $(LIBDIR)/*.so $(FILTERDIR)/*~ $(FILTERDIR)/*.o

distclean: clean
	

install:
	@printf "Installing MeteoIO\n"

documentation:
	doxygen $(SRCDIR)/config.dox

####### Compile

.cc.o: $*.cc $*.h
	$(CXX) $(CCFLAGS) -c $< $(INCLUDE) -o $@

$(LIBDIR)/libmeteoIO.a: $(METEOIO_OBJ)
	ar -r $@ $(METEOIO_OBJ)
	ranlib $@

$(LIBDIR)/libfilter.a: $(FILTER_OBJ)
	ar -r $@ $(FILTER_OBJ)
	ranlib $@

$(LIBDIR)/libBoschungIO.so: $(SRCDIR)/BoschungIO.cc $(SRCDIR)/BoschungIO.h $(LIBDIR)/libmeteoIO.a $(LIBDIR)/libfilter.a
ifeq ($(BOSCHUNGIO),yes)
	$(CXX) $(CCFLAGS) -fPIC $(INCLUDE) $(shell pkg-config --cflags libxml++-2.6) -c -o $(SRCDIR)/BoschungIO.o $(SRCDIR)/BoschungIO.cc 
	$(CXX) $(CCFLAGS) -rdynamic -shared -Wl,-soname,libBoschungIO.so -o $@ $(SRCDIR)/BoschungIO.o \
	$(LDFLAGS_SEQ) $(LDFLAGS) $(shell pkg-config --libs libxml++-2.6)
endif

$(LIBDIR)/libImisIO.so: $(SRCDIR)/ImisIO.cc $(SRCDIR)/ImisIO.h $(LIBDIR)/libmeteoIO.a $(LIBDIR)/libfilter.a
ifeq ($(IMISIO),yes)
	g++-3.3 -Wall $(DEBUG) -fPIC $(INCLUDE) $(INCLUDE_ORA) $(SRCDIR)/ImisIO.cc -c -o $(SRCDIR)/ImisIO.o
	g++-3.3 $(DEBUG) -rdynamic -shared -Wl,-soname,libImisIO.so -o $@ $(SRCDIR)/ImisIO.o \
	$(LDFLAGS_SEQ) $(LDFLAGS) $(LIBS_ORA) $(LDFLAGS_ORA)
endif

meteoIO.module: libmeteoOIOparoc.a $(SRCDIR)/PackMeteoIO_par.o
	$(PAROCC) $(CCFLAGS) -object -parocld=$(LINKER) -o $@ $(METEOIODIR)/PackMeteoIO_par.o $(LDFLAGS) $(LDFLAGS_PAROC)

libmeteoIOparoc.a:  $(METEOIO_OBJ_PAROC)
	ar -r $@ $(METEOIO_OBJ_PAROC)
	ranlib $@

libfilterparoc.a: $(FILTER_OBJ_PAROC)
	ar -r $@ $(FILTER_OBJ_PAROC)
	ranlib $@

