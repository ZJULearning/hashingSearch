GXX=g++ -std=c++11
OPTM=-O3 -march=native -fopenmp
CPFLAGS=$(OPTM) -Wall 
LDFLAGS=$(OPTM) -Wall -lboost_timer -lboost_chrono -lboost_system 

INCLUDES=-I./ -I./algorithm -I./general

SAMPLES=$(patsubst %.cc, %, $(wildcard samples/*.cc samples_hashing/*.cc))
SAMPLE_OBJS=$(foreach sample, $(SAMPLES), $(sample).o)

HEADERS=$(wildcard ./*.hpp ./*/*.hpp)

all: $(SHARED_LIB) $(SAMPLES)


$(SAMPLES): %: %.o
	$(GXX) $^ -o $@ $(LDFLAGS) $(LIBS)

%.o: %.cpp $(HEADERS)
	$(GXX) $(CPFLAGS) $(INCLUDES) -c $*.cpp -o $@

%.o: %.cc $(HEADERS)
	$(GXX) $(CPFLAGS) $(INCLUDES) -c $*.cc -o $@

clean:
	rm -rf $(OBJS)
	rm -rf $(SHARED_LIB)
	rm -rf $(SAMPLES)
	rm -rf $(SAMPLE_OBJS)
