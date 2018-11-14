CC = g++
CFLAGS = -g -funroll-loops -std=c++17 -O2 -Wall -pedantic -llapack
FILE = qp
SRCS = $(FILE).h $(FILE).cpp
OBJS = $(FILE).o
TEST = test
TESTSRC = $(TEST).cpp

all: product test

product: $(SRCS) \
; $(CC) $(CFLAGS) -c $(SRCS)

test: $(TESTSRC) \
; $(CC) $(CFLAGS) -o $(TEST) $(TESTSRC) $(OBJS)

clean: \
; $(RM) $(OBJS) $(TEST)

write: \
; echo $(OBJS)
