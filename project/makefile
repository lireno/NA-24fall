# Compiler and flags
CXX = g++
CXXFLAGS =-std=c++17
HEADERS = $(wildcard ./src/*.h)

# Directories
SRC_DIR = ./src
BIN_DIR = ./bin
TEST_DIR = ./test
PLOT_DIR = ./plotting
DOC_DIR = ./doc

# Files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(BIN_DIR)/%.o, $(SOURCES))
TEST_SOURCES = $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJECTS = $(patsubst $(TEST_DIR)/%.cpp, $(BIN_DIR)/%.o, $(TEST_SOURCES))

.PHONY: all run report test plot clean

# Default target
all: run report

# Compile and run the source files
run: $(OBJECTS)
	@echo "Compiling and running source files..."
	$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/main $(OBJECTS)
	@$(BIN_DIR)/main

$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Generate reports
report: $(DOC_DIR)/design.pdf $(DOC_DIR)/report.pdf

$(DOC_DIR)/%.pdf: $(DOC_DIR)/%.tex
	@echo "Generating $@..."
	pdflatex -output-directory=$(DOC_DIR) $<

# Compile and run tests
test: ./test/test.cpp $(HEADERS)
	@echo "Compiling and running tests..."
	g++ -g -O0 -I./src -o t test/testcurve.cpp
	gdb t

$(BIN_DIR)/%.o: $(TEST_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Generate plots
plot: $(PLOT_DIR)/plot.cpp  $(HEADERS)
	@echo "Generating plots..."
	$(CXX) $(CXXFLAGS) -I./src -o $(PLOT_DIR)/plot $(PLOT_DIR)/plot.cpp
	$(PLOT_DIR)/plot
	rm $(PLOT_DIR)/plot
	python3 $(PLOT_DIR)/plot.py


# Clean up generated files
clean:
	@echo "Cleaning up generated files..."
	rm -rf $(BIN_DIR)/*.o $(BIN_DIR)/main $(BIN_DIR)/test
	rm -f $(DOC_DIR)/*.pdf
	rm -rf $(DOC_DIR)/*.aux $(DOC_DIR)/*.log
	rm -f $(PLOT_DIR)/plot
	rm -f $(PLOT_DIR)/data/*.txt $(PLOT_DIR)/figures/*.png

