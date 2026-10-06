CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -fopenmp
TARGET := viper_ct
SOURCE := src/viper_ct.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) $(SOURCE) -o $(TARGET)

clean:
	rm -f $(TARGET)
