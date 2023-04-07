CC = g++
CFLAGS = -g -std=c++17 -O2 -Wall -pedantic -lblas -llapack
FILE = QP
SRCS = $(FILE).h $(FILE).cpp
OBJS = $(FILE).o MathUtil.o LinearSolver.o
MATH = MathUtil.h MathUtil.cpp
LAPACK = LinearSolver.h LinearSolver.cpp
TEST = test
TESTSRC = $(TEST).cpp

all: lapack product test

product: $(SRCS) \
; $(CC) -c $(MATH) $(SRCS) $(CFLAGS)

math: $(MATH) \
; $(CC) -c $(MATH) $(CFLAGS)

lapack: $(LAPACK) \
; $(CC) -c $(LAPACK) $(CFLAGS)

test: $(TESTSRC) \
; $(CC) -o $(TEST) $(TESTSRC) $(OBJS) $(CFLAGS)

clean: \
; $(RM) $(OBJS) $(TEST)

write: \
; echo $(OBJS)
