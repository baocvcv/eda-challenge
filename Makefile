CC := g++
FLAGS := -std=c++11
OPTFLAGS := -O3
DBG_FLAGS := -std=c++11 -g -DDebug
source := main.cpp
headers := $(wildcard *.h *.hpp)

.PHONY: all
all: main

main: $(source) $(headers)
	$(CC) $(FLAGS) $(OPTFLAGS) $< -o $@

check: $(source) $(headers)
	$(CC) $(FLAGS) -O0 $< -o $@
	rm -rf $@

debug: $(source) $(headers)
	$(CC) $(DBG_FLAGS) $< -o $@

.PHONY: run
run: main
	./main ../data/testdata-1

.PHONY: clean
clean:
	rm -rf main main.exe debug debug.exe