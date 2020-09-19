LD_LIBRARY_PATH := /home/esash/Documents/Projects/CellularAutomata/CellularAutomata_Fast/lib

out: libcaf.so main.o
	if [ ! -d "./bin" ]; then mkdir bin; fi
	g++ -g -o bin/out obj/main.o obj/plot.o -L$(LD_LIBRARY_PATH) -lcaf -ldl

main.o:
	if [ ! -d "./obj" ]; then mkdir obj; fi
	g++ -g -o obj/main.o -c src/main.cpp -Iinclude/
	g++ -g -o obj/plot.o -c src/plot.cpp -Iinclude/

libcaf.so: 
	g++ -g -o lib/libcaf.so -Wall -fPIC -Ilib/include/ -Ilib/external -shared \
	lib/src/log.cpp lib/src/network.cpp \
	lib/src/population.cpp lib/src/util.cpp lib/src/config.cpp \
	lib/src/spike_generator.cpp

clean:
	rm -rf bin obj lib/libcaf.so results/*

test:
	make testSpikeGenerator
	make testConfig
	make testPopulation
	make testModel
	./bin/test/testSpikeGenerator
	./bin/test/testConfig
	./bin/test/testPopulation
	./bin/test/testModel

testSpikeGenerator:
	if [ ! -d "./bin" ]; then mkdir bin; fi
	if [ ! -d "./bin/test" ]; then mkdir bin/test; fi
	g++ -o bin/test/testSpikeGenerator lib/src/spike_generator.cpp lib/src/util.cpp \
	lib/test/TestSpike_Generator.cpp -lcppunit

testConfig:
	if [ ! -d "./bin" ]; then mkdir bin; fi
	if [ ! -d "./bin/test" ]; then mkdir bin/test; fi
	g++ -o bin/test/testConfig lib/src/config.cpp lib/src/util.cpp \
	lib/src/log.cpp lib/src/spike_generator.cpp \
	lib/src/network.cpp \
	lib/src/population.cpp lib/test/TestConfig.cpp -lcppunit

testPopulation:
	if [ ! -d "./bin" ]; then mkdir bin; fi
	if [ ! -d "./bin/test" ]; then mkdir bin/test; fi
	g++ -o bin/test/testPopulation lib/test/TestPopulation.cpp -lcppunit -L$(LD_LIBRARY_PATH) -lcaf -ldl

testModel:
	if [ ! -d "./bin" ]; then mkdir bin; fi
	if [ ! -d "./bin/test" ]; then mkdir bin/test; fi
	g++ -o bin/test/testModel lib/test/TestModel.cpp -lcppunit -L$(LD_LIBRARY_PATH) -lcaf -ldl
