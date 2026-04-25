# ================================================================
#  Makefile for C++ MPI + OpenMP Project
#  Automatically compiles all .cpp files in src/
#  Automatically tracks all headers in include/
#  Creates build/ and bin/ directories
#  Supports: make, make run, make clean
# ================================================================

# --- Compiler Settings ---
CXX     := mpicxx
CXXFLAGS := -O2 -Wall -std=c++23 -fopenmp

# --- Include & Library Paths ---
INCLUDES := -Iinclude
LIBS     := -lm

# --- Directories ---
SRC_DIR  := src
INC_DIR  := include
BUILD_DIR := build
BIN_DIR  := bin

# --- Output Executable ---
TARGET := $(BIN_DIR)/Program

# --- Automatically find all source and header files ---
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
HEADERS := $(wildcard $(INC_DIR)/*.hpp) $(wildcard $(INC_DIR)/*.h)

# --- Object files go into build/ ---
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))

# --- Parallel run defaults ---
NP  ?= 4        # default MPI ranks
OMP ?= 1        # default OpenMP threads

# ================================================================
#  Rules
# ================================================================

# Default rule builds the main executable
all: $(TARGET)

# Rule to build final executable
$(TARGET): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
	@echo "Built executable: $(TARGET)"

# Rule to build each object file
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# ================================================================
#  Run the MPI program
# ================================================================
run: $(TARGET)
	OMP_NUM_THREADS=$(OMP) mpirun -np $(NP) $(TARGET) $(ARGS)

# ================================================================
#  Clean
# ================================================================
clean:
	rm -rf $(BUILD_DIR)/*.o $(BIN_DIR)/Program
	@echo "Cleaned build files."
