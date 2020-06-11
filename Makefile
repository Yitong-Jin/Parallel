MAKE = MAKE
TARGET = my_code
SOURCE = parallel.cpp

default:
	mpicxx -o $(TARGET) $(SOURCE) --std=c++11