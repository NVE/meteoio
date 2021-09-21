CC      = g++
CFLAGS  = -Wall -Wextra
DEBUG   = -g
METEOIODIR = /usr/local/include

LIBS    = -rdynamic -lstdc++ -lmeteoio -ldl
INCLUDE = -I. -I$(METEOIODIR)

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
all: grid_interpolations

grid_interpolations: grid_interpolations.o
	$(CC) $(DEBUG) -o $@ $@.o ${LIBS}

clean:
	rm -rf *~ *.o grid_interpolations

