OBJS = compareHistograms.o
EXE = compareHistograms

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTLIBS     += -lEG

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++ -W -Wall -Wextra -pedantic -std=c++11
FLAGS = -Wall -g $(INCFLAGS) $(LDFLAGS) -fPIC

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE) $(EXE2)

$(EXE): $(OBJS)
	$(CXX) -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS)

$(EXE2): $(OBJS2)
	$(CXX) -o $(EXE2) $(OBJS2) $(ROOTFLAGS) $(ROOTLIBS)

%.o: %.cxx
	$(COMPILE)  $< 

clean:
	rm -f *.o *.so compareHistograms
