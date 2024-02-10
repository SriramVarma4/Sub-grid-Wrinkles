CC = g++
CFLAGS = -Wall -std=c++17 -I ./eigen-3.4.0
LDFLAGS = -lGL -lGLU -lglut -lfreeimage -w

SRC = main.cpp
SRC2 = main2.cpp
OUT = video
OUT2 = video2

all: $(OUT) $(OUT2)

$(OUT): $(SRC)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

$(OUT2): $(SRC2)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

run: $(OUT)
	./$(OUT)

run2: $(OUT2)
	./$(OUT2)

clean:
	rm -f $(OUT) $(OUT2)
