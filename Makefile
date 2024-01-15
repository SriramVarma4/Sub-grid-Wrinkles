CC = g++
CFLAGS = -Wall -std=c++11 -w
LDFLAGS = -lGL -lGLU -lglut -lfreeimage

SRC = main.cpp
OUT = video

all: $(OUT)

$(OUT): $(SRC)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

run: $(OUT)
	./$(OUT)

clean:
	rm -f $(OUT)
