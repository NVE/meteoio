#
# Makefile for the METEOIO Library
#
#################################################
# USER OPTIONS: MAKE CHANGES IN THIS SECTION ONLY
#################################################
####### USER CONFIGURATION
#Destination system: either zeus, grid or safe
DEST	 = safe
#DEST	 = grid
#DEST	 = zeus

# build mode: release or debug
MODE	 = release
MODE	 = debug

# plugins to build (yes/no)
BOSCHUNGIO 	= no
IMISIO		= no

####### COMPILERS AND OPTIONS
#SVNREV_GEN	= $(shell main/version.sh)
#SVNREV		= $(eval SVNREV := $(SVNREV_GEN))$(SVNREV)

ifeq ($(DEST),safe)		#safe defaults
	CC       = gcc -DGNU
	CXX      = g++
	FF	 = gfortran
	LINKER   = g++ -DGNU
	PAROCC   = parocc
	ARCH     = -march=pentium3
	OPTIM    = -fomit-frame-pointer -O3
	FFOPTIM  = $(OPTIM) -march=pentium3
	FFEXTRA  = $(FFOPTIM)
	#MAKE_OPTIM = -pipe
	DEBUG    = -g -O0 -D__DEBUG #-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC #-dr 
	WARNING  = -Wall #-Weffc++ #-pedantic-errors #-Werror
endif
ifeq ($(DEST),grid)		#for grid
	CC	 = colorgcc
	CXX	 = colorgcc
	FF	 = gfortran
	LINKER   = g++ -DGNU
	PAROCC   = parocc
	ARCH     = -march=native #-mmmx -msse -msse2 -mfpmath=sse -malign-double
	OPTIM    = -fomit-frame-pointer -O3 #-fdata-sections
	FFOPTIM  = $(OPTIM) -march=native #-m3dnow -m64 -mmmx -msse -msse2 -ffast-stdlib 
	FFEXTRA  = $(FFOPTIM) #-fno-second-underscore
	MAKE_OPTIM = -combine -pipe
	DEBUG    = -g -O0 -D__DEBUG #-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC #-dr 
	WARNING  = -Wall -Wextra #-Weffc++ #-pedantic-errors #-Werror
endif
ifeq ($(DEST),zeus)		#for Zeus
	CC       = gcc -DGNU
        #CC      = pathcc -DGNU
        CXX      = g++
        #CXX     = pathcc -DGNU
        FF       = pathf90 -DGNU
        LINKER   = g++
        PAROCC   = parocc
        ARCH     = -march=x86-64 #-mmmx -msse -msse2 -m3dnow -mfpmath=sse #-malign-double
        #ARCH     = -march=opteron -mmmx -msse -msse2 -m3dnow
        OPTIM    = #-fomit-frame-pointer #-O3 -fdata-sections
        #OPTIM   = -O2 -Wstrict-aliasing #-fno-strict-aliasing
        FFOPTIM  = $(OPTIM) -march=opteron #-m3dnow -m64 -mmmx -msse -msse2 -ffast-stdlib 
        FFEXTRA  = -fno-second-underscore $(FFOPTIM)
	MAKE_OPTIM = -pipe
	DEBUG    = -g -O0 -D__DEBUG
	WARNING  = -Wall
endif

ifeq ($(MODE),debug)		#for debug
	CFLAGS   = $(MAKE_OPTIM) $(ARCH) $(DEBUG) $(WARNING)
	CCFLAGS  = $(MAKE_OPTIM) $(ARCH) $(DEBUG) $(WARNING)
	FFLAGS	 = $(MAKE_OPTIM) $(ARCH) $(DEBUG) $(WARNING)
else				#for release
	CFLAGS   = $(MAKE_OPTIM) $(ARCH) $(OPTIM) $(WARNING)
	CCFLAGS  = $(MAKE_OPTIM) $(ARCH) $(OPTIM) $(WARNING)
	FFLAGS	 = $(MAKE_OPTIM) $(ARCH) $(FFOPTIM) $(WARNING)
endif
#CCFLAGS += -DSVNREV="\"$(SVNREV)\""

####### DIRECTORIES
TARGET		= meteoio
TARGET_PAROC	= meteoio_paroc

SRCDIR		= ./src
LIBDIR		= ./lib
FILTERDIR	= $(SRCDIR)/filter
TOOLSDIR	= ./tools

#################################################
# END USER OPTIONS
#################################################

LIBS		= -lc -ldl
LDFLAGS_SEQ	= -L$(LIBDIR) -lmeteoio -lfilter
LD_PAROC	= -L$(LIBDIR) -lmeteoioparoc -lfilterparoc
LDFLAGS_PAROC	= $(LD_PAROC)
LDFLAGS		= $(LIBS) -rdynamic
INCLUDE		= -I$(SRCDIR) -I$(FILTERDIR)

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
		$(SRCDIR)/Laws.o \
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
		$(SRCDIR)/IOHandler.stub.o \
		$(SRCDIR)/IOHandler_par.o \
		$(SRCDIR)/A3DIO_par.o \
		$(SRCDIR)/LegacyIO.stub.o \
		$(SRCDIR)/MeteoData_par.o \
		$(SRCDIR)/StationData_par.o \
		$(SRCDIR)/ConfigReader_par.o \
		$(SRCDIR)/Date_IO_par.o \
		$(SRCDIR)/Grid2DObject_par.o \
		$(SRCDIR)/IOExceptions_par.o \
		$(SRCDIR)/IOUtils_par.o \
		$(SRCDIR)/Laws_par.o \
		$(SRCDIR)/libinterpol2D_par.o \
		$(SRCDIR)/libinterpol1D_par.o \
		$(SRCDIR)/Meteo1DResampler_par.o \
		$(SRCDIR)/Meteo2DInterpolator_par.o \
		$(SRCDIR)/MeteoBuffer_par.o \
		$(SRCDIR)/Meteo1DResampler_par.o \
		$(SRCDIR)/DynamicLibrary_par.o \
		$(SRCDIR)/marshal_meteoio_par.o 

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

TOOLS_OBJ =	$(TOOLSDIR)/createA3DFiles.o

####### Build rules SEQ

all: seq 
	@printf "*** MeteoIO compiled for \033[36m%s\033[0m\n" $(DEST)

help:
	@printf "MeteoIO Makefile targets for \033[36m%s\033[0m:\n" $(DEST)
	@printf " \033[36mseq\033[0m \n"
	@printf " \033[36mcreateA3DFiles\033[0m \n"
	@printf " \033[36mparoc\033[0m \n"
	@printf " \033[36minstall\033[0m \n"
	@printf " \033[36mclean\033[0m \n"
	@printf " \033[36mdistclean\033[0m \n"
	@printf " \033[36mdocumentation\033[0m\n"

seq: $(LIBDIR)/libmeteoio.a build_dynamiclibs build_staticlibs

build_staticlibs: $(LIBDIR)/libfilter.a

build_dynamiclibs: $(LIBDIR)/libboschungio.so $(LIBDIR)/libimisio.so

############## PAROC ##############
paroc: meteoIO_lib_paroc meteoIO_module_paroc filter_lib_paroc

filter_lib_paroc: $(LIBDIR)/libfilterparoc.a

meteoIO_lib_paroc: $(LIBDIR)/libmeteoioparoc.a 

meteoIO_module_paroc: $(LIBDIR)/meteoio.module
##############  END  ##############

clean:
	rm -f $(SRCDIR)/*~ $(SRCDIR)/*.o $(FILTERDIR)/*~ $(FILTERDIR)/*.o $(TOOLSDIR)/*.o

distclean: clean
	rm $(TOOLSDIR)/createA3DFiles
	rm $(LIBDIR)/*.a $(LIBDIR)/*.so $(LIBDIR)/*.module

install:
	@printf "**** Installing MeteoIO\n"

documentation:
	doxygen $(SRCDIR)/config.dox

createA3DFiles: $(TOOLS_OBJ)
	$(CXX) $(CCFLAGS) $< $(INCLUDE) $(LDFLAGS_SEQ) $(LDFLAGS) -o $(TOOLSDIR)/$@

####### Compile

.cc.o: $*.cc $*.h
	$(CXX) $(CCFLAGS) -c $< $(INCLUDE) -o $@

%_par.o : %.cc
	$(PAROCC) $(CCFLAGS) -c $< $(INCLUDE) -o $@

%.stub.o : %.ph
	$(PAROCC) $(CCFLAGS) -c $< $(INCLUDE) -o $@


$(LIBDIR)/libmeteoio.a: $(METEOIO_OBJ)
	@printf "**** Compiling libmeteoIO\n"
	ar -r $@ $(METEOIO_OBJ)
	ranlib $@

$(LIBDIR)/libfilter.a: $(FILTER_OBJ)
	ar -r $@ $(FILTER_OBJ)
	ranlib $@

$(LIBDIR)/libboschungio.so: $(SRCDIR)/BoschungIO.cc $(SRCDIR)/BoschungIO.h $(LIBDIR)/libmeteoio.a $(LIBDIR)/libfilter.a
ifeq ($(BOSCHUNGIO),yes)
	@printf "**** Compiling Boschung plugin\n"
	$(CXX) $(CCFLAGS) -fPIC $(INCLUDE) $(shell pkg-config --cflags libxml++-2.6) -c -o $(SRCDIR)/BoschungIO.o $(SRCDIR)/BoschungIO.cc 
	$(CXX) $(CCFLAGS) -rdynamic -shared -Wl,-soname,libboschungio.so -o $@ $(SRCDIR)/BoschungIO.o $(LDFLAGS_SEQ) $(LDFLAGS) $(shell pkg-config --libs libxml++-2.6)
endif

$(LIBDIR)/libimisio.so: $(SRCDIR)/ImisIO.cc $(SRCDIR)/ImisIO.h $(LIBDIR)/libmeteoio.a $(LIBDIR)/libfilter.a
ifeq ($(IMISIO),yes)
	@printf "**** Compiling Imis plugin\n"
	$(CXX) $(CCFLAGS) -fPIC $(INCLUDE) -I$(ORACLE_HOME)/rdbms/public $(SRCDIR)/ImisIO.cc -c -o $(SRCDIR)/ImisIO.o
	$(CXX) $(CCFLAGS) -rdynamic -shared -Wl,-rpath,$(ORACLE_HOME)/lib,-soname,libimisio.so -o $@ $(SRCDIR)/ImisIO.o $(LDFLAGS_SEQ) $(LDFLAGS) -L$(ORACLE_HOME)/lib -locci -lclntsh -lstdc++
endif

$(LIBDIR)/meteoio.module: $(LIBDIR)/libmeteoioparoc.a $(LIBDIR)/libfilterparoc.a $(SRCDIR)/PackMeteoIO_par.o
	$(PAROCC) $(CCFLAGS) -object -parocld=$(LINKER) -o $@ $(SRCDIR)/PackMeteoIO_par.o $(LDFLAGS) $(LDFLAGS_PAROC)

$(LIBDIR)/libmeteoioparoc.a:  $(METEOIO_OBJ_PAROC)
	ar -r $@ $(METEOIO_OBJ_PAROC)
	ranlib $@

$(LIBDIR)/libfilterparoc.a: $(FILTER_OBJ_PAROC)
	ar -r $@ $(FILTER_OBJ_PAROC)
	ranlib $@

