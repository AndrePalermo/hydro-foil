OBJ_DIR = obj
SRC_DIR = src

CXX = g++
CXXFLAGS = -fPIC -O3 -std=c++20 -march=native -fopenmp -lm
OPEN_MP_FLAG = -DOPEN_MP #remove if you don't want to use openmp

SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))

TARGET = foil #Freeze-Out IntegraL

########################################################################################

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OPEN_MP_FLAG)$^ -o $@
	@echo "Compilation successful!"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(OPEN_MP_FLAG) -c $< -o $@

