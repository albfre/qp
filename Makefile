CC = g++
CFLAGS = -g -std=c++17 -O2 -Wall -pedantic -llapack
FILE = qp
SRCS = $(FILE).h $(FILE).cpp
OBJS = $(FILE).o MathUtil.o
MATH = MathUtil.h MathUtil.cpp
TEST = test
TESTSRC = $(TEST).cpp

all: product test

product: $(SRCS) \
; $(CC) $(CFLAGS) -c $(MATH) $(SRCS)

math: $(MATH) \
; $(CC) $(CFLAGS) -c $(MATH)

test: $(TESTSRC) \
; $(CC) $(CFLAGS) -o $(TEST) $(TESTSRC) $(OBJS)

clean: \
; $(RM) $(OBJS) $(TEST)

write: \
; echo $(OBJS)
