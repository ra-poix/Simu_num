CC=g++
CFLAGS=-Wall -std=c++11
EXEC=exec
SRC= $(wildcard Option/*.cpp Model/*.cpp Generateur/*.cpp  *.cpp)
OBJ= $(SRC:Option/*.cpp=object/*.o Model/*.cpp=object/*.o Generateur/*.cpp=object/*.o .cpp=.o )

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY: clean mrproper

clean:
	@rm -rf *.o

mrproper: clean
	@rm -rf $(EXEC)