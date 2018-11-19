CC = g++
CFLAGS = -g -std=c++17 -O2 -Wall -pedantic -llapack
FILE = QP
SRCS = $(FILE).h $(FILE).cpp
OBJS = $(FILE).o MathUtil.o LinearSolver.o
MATH = MathUtil.h MathUtil.cpp
LAPACK = LinearSolver.h LinearSolver.cpp
TEST = test
TESTSRC = $(TEST).cpp

all: lapack product test

product: $(SRCS) \
; $(CC) $(CFLAGS) -c $(MATH) $(SRCS)

math: $(MATH) \
; $(CC) $(CFLAGS) -c $(MATH)

lapack: $(LAPACK) \
; $(CC) $(CFLAGS) -c $(LAPACK)

test: $(TESTSRC) \
; $(CC) $(CFLAGS) -o $(TEST) $(TESTSRC) $(OBJS)

clean: \
; $(RM) $(OBJS) $(TEST)

write: \
; echo $(OBJS)
