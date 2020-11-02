CC := g++
FLAGS := -std=c++14 -O3
DBG_FLAGS := -std=c++14 -g -DDebug
source := main.cpp
headers := $(wildcard *.h *.hpp)

main: $(source) $(headers)
	$(CC) $(FLAGS) $< -o $@

debug: $(source) $(headers)
	$(CC) $(DBG_FLAGS) $< -o $@

.PHONY: run
run: main
	./main ../data/testdata-1