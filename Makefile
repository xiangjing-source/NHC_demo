CXX = g++
CXXFLAGS = -O2 -std=c++11

all: md

md: md.cpp
	@$(CXX) $(CXXFLAGS) md.cpp -o md
	@echo "Compilation successful."

clean:
	@rm -f md *.o
	@echo "Cleaned build files."
