CC=gcc
CFLAGS=-Wall -fopenmp -std=c99 -O3 -g -I.
LDFLAGS=-lm
SOURCES=$(wildcard *.c)
OBJECTS=$(SOURCES:.c=.o)
TARGET=distr_5
	
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
	
clean:
	rm -r -f *.o $(TARGET)
