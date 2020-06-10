# the compiler: gcc for C program, define as g++ for C++
CC = gcc
#  compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall -lpthread -lm

# the build target executable:
TARGET = srf

all: main.c
	$(CC) *.c ./lib/* $(CFLAGS) -o $(TARGET) 

clean:
	$(RM) $(TARGET)
