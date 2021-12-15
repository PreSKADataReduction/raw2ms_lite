target=bin/raw2ms_splited

all:$(target)

INC=-I ../mscreate/include -I /usr/include/casacore/ -I ./include
CXXFLAGS=-O3 -std=c++11
LDFLAGS=-lcasa_casa -lcasa_ms -lcasa_measures -lcasa_tables

bin/raw2ms_splited:obj/raw2ms_splited.o obj/date_time.o obj/mscreate.o
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/raw2ms_splited.o:src/raw2ms_splited.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@


obj/date_time.o:src/date_time.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

obj/mscreate.o: src/mscreate.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

install: $(target)
	cp -r bin/* /usr/local/bin/

clean:
	rm -f `find -iname *.o` `find -iname *~`
