objects = objects/laguerre.o objects/big_numbers.o objects/precomputed_laguerre_roots.o objects/PyLaplace.o
options = -O2 -fPIC -std=c++17 -pipe
INSTALL = /usr/local/lib/python3.5/dist-packages/
PY_INCLUDE = /usr/include/python3.5/
CC = g++-7
lib/laplace.so: $(objects) 
	$(CC) -O2 -fPIC -shared -o lib/laplace.so $(objects)
objects/big_numbers.o: big_numbers.h big_numbers.cpp
	$(CC) $(options) -c big_numbers.cpp -o objects/big_numbers.o

objects/laguerre.o: objects/big_numbers.o polynomials.h laguerre.h laguerre.cpp
	$(CC) $(options) -c laguerre.cpp -o objects/laguerre.o

objects/precomputed_laguerre_roots.o: precomputed_laguerre_roots.h objects/big_numbers.o precomputed_laguerre_roots.cpp
	$(CC) $(options) -c precomputed_laguerre_roots.cpp -o objects/precomputed_laguerre_roots.o

objects/PyLaplace.o: laplace.h PyLaplace.cpp
	$(CC) $(options) -c -I$(PY_INCLUDE) PyLaplace.cpp -o objects/PyLaplace.o
clean:
	rm $(objects) lib/laplace.so
install: $(lib/laplace.so)
	cp lib/laplace.so $(INSTALL)
