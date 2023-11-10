target=bin/raw2ms_splited bin/ms2uvmap

all:$(target)

INC=-I /usr/include/casacore/ -I ./include
CXXFLAGS=-O3 -std=c++11
LDFLAGS=-lcasa_casa -lcasa_ms -lcasa_measures -lcasa_tables -lcfitsio -lboost_program_options -lboost_system

OBJS=obj/date_time.o obj/fio.o obj/fitsfile.o obj/fits_trait.o obj/region.o obj/region_imp.o obj/mscreate.o

bin/raw2ms_splited:obj/raw2ms_splited.o $(OBJS)
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g


bin/ms2uvmap: obj/ms2uvmap.o $(OBJS)
	mkdir -p bin&&$(CXX) -o $@ $^ $(LDFLAGS) -g

obj/raw2ms_splited.o:src/raw2ms_splited.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

obj/ms2uvmap.o:src/ms2uvmap.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

obj/fio.o: src/fio.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

obj/fitsfile.o: src/fitsfile.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

obj/fits_trait.o: src/fits_trait.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

obj/region.o: src/region.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

obj/region_imp.o: src/region_imp.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

obj/date_time.o:src/date_time.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

obj/mscreate.o: src/mscreate.cpp
	mkdir -p obj&&$(CXX) -c $< $(CXXFLAGS) $(INC) -g -o $@

install: $(target)
	cp -r bin/* /usr/local/bin/

clean:
	rm -f `find -iname *.o` `find -iname *~`
