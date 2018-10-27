EXT = glx64
MEX = /usr/local/MATLAB/R2016a/bin/mex

BASE = ${shell dirname `pwd` }
CFLAGS = -Wall 
MEXFLAGS = -lmex  -lm

all: ../triangle.mex${EXT} ../refine.mex${EXT} ../triangle
test: test_triangle

../triangle: tristruct.o util.o trilib.c
	${CC} ${CFLAGS} $^ -o $@ -lm

../triangle.mex${EXT}: triangle.c tristruct.o util.o trilib.o
	${MEX} ${MEXFLAGS} CFLAGS="${CFLAGS}" $^ 
	mv triangle.mex${EXT} ../triangle.mex${EXT}

../refine.mex${EXT}: refine.c tristruct.o util.o trilib.o
	echo '#define BASE "${BASE}"' > base.h
	${MEX} CFLAGS="${CFLAGS}" $^ 
	mv refine.mexmaci64 ../refine.mexmaci64

test_triangle: test_triangle.o tristruct.o util.o trilib.o
	${CC} ${CFLAGS} test_triangle.o tristruct.o util.o trilib.o  -lm -o $@

trilib.o: trilib.c
	${CC} ${CFLAGS} -c -g -DTRILIBRARY $<

%.o: %.cpp
	${CC} ${CFLAGS} -c -g $<

clean:
	\rm -f  ../triangle.mex${EXT} ../refine.mex${EXT} test_triangle ../triangle *.o *~ *.out

# dependances

triangle.o: triangle.c triangle.h tristruct.h util.h
trilib.o: trilib.c triangle.h
tristruct.o: tristruct.c tristruct.h triangle.h util.h
util.o: util.c util.h

