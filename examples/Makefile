ALGLIB=alglib
LALGLIB=-L$(ALGLIB) -lalglib
IALGLIB=-I$(ALGLIB)

CC=g++
CFLAGS=-I../src -Itk -Ilibspline -I${ALGLIB} -Wall -O2 -std=c++11
LIB=libspline

all: basics bench lib

basics: basics.cpp
	${CC} ${CFLAGS} -c $< -o basics.o
	$(CC) $(CFLAGS) basics.o -o $@ $(LDFLAGS)

bench: bench.cpp
	${CC} ${CFLAGS} -c $< -o bench.o
	${CC} $(CFLAGS) bench.o -o $@ ${LALGLIB}

lib: libspline/$(LIB).cpp
	$(CC) $(CFLAGS) -fPIC -Wextra -c $^ -o libspline/$(LIB).o
	$(CC) -shared -o libspline/$(LIB).so libspline/$(LIB).o

%.o: %.cpp
	${CC} ${CFLAGS} -c $< -o $@

clean:
	-rm -f basics bench *.o *.so libspline/*.o libspline/*.so 2>/dev/null