#!/bin/bash

gcc -I./src -M ./examples/maths/*.cpp ./examples/finance/*.cpp \
	./examples/physics/*.cpp ./tests/*.cpp \
       	| sed -e 's/[\\ ]/\n/g' | \
        sed -e '/^$/d' -e '/\.o:[ \t]*$/d' | \
        ctags -L - --c++-kinds=+p --fields=+iaS --extra=+f+q
