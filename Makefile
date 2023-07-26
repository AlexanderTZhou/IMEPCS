CC = g++
CFLAGS = -g -Wall -std=c++11
TARGET = epcli

all: $(TARGET)

$(TARGET) : $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp
