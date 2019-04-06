CC=gcc
CFLAGS=-Wall
EXEC=exec
OBJ=main

all: $(EXEC)

%.o: %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

$(EXEC): main.o Call.o
	$(CC) -o $@ $< 

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)