all: compile

compile: src/tiny_ekf.c src/vecmat_utils.c
	gcc -fPIC -I./inc -c -fPIC src/tiny_ekf.c -o build/tinyekf.o -lm
	gcc -fPIC -I./inc -c -fPIC src/vecmat_utils.c -o build/vecmat_utils.o -lm
	ar -rcs build/libtinyekf.a build/tinyekf.o build/vecmat_utils.o

clean:
	rm -rf build/*

