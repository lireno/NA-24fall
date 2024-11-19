# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -I./src

# Directories
SRC_DIR = ./src
BIN_DIR = ./bin
TEST_DIR = ./test
PLOT_DIR = ./plotting
DOC_DIR = ./doc

# Files
HEADERS = $(wildcard $(SRC_DIR)/*.h)
SOURCES = $(wildcard assignment/*.cpp)
OBJECTS = $(patsubst assignment/%.cpp, $(BIN_DIR)/%, $(SOURCES))
TEST_SOURCES = $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJECTS = $(patsubst $(TEST_DIR)/%.cpp, $(BIN_DIR)/%, $(TEST_SOURCES))
PLOT_SOURCES = $(PLOT_DIR)/plot.cpp

# Targets
.PHONY: all run test report data plot plot1 clean

# Default target
all: $(OBJECTS) report

# Build assignment executables
$(BIN_DIR)/%: assignment/%.cpp $(HEADERS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN_DIR)/%: $(BIN_DIR)/%
	$(CXX) $(CXXFLAGS) -o $@ $<
	@echo "Built $@"

run: $(OBJECTS)
	@for prog in $(OBJECTS); do \
		exec_file=$(BIN_DIR)/$$(basename $$prog); \
		echo "Running $$exec_file..."; \
		$$exec_file; \
		echo ""; \
	done

# Build and run tests
test: $(TEST_SOURCES) $(HEADERS)
	@echo "Compiling and running tests..."
	$(CXX) $(CXXFLAGS) -o $(BIN_DIR)/tests $(TEST_SOURCES)
	$(BIN_DIR)/tests

# Generate data from plot.cpp
data: $(PLOT_SOURCES)
	@echo "Generating data..."
	$(CXX) $(CXXFLAGS) -o $(PLOT_DIR)/plot $(PLOT_SOURCES)
	$(PLOT_DIR)/plot
	rm -f $(PLOT_DIR)/plot

# Generate plots using Python scripts
plot:
	python3 $(PLOT_DIR)/plot.py

plot1:
	python3 $(PLOT_DIR)/plot_in_one_figure.py

# Build LaTeX reports
report: $(DOC_DIR)/design.pdf $(DOC_DIR)/report.pdf

$(DOC_DIR)/%.pdf: $(DOC_DIR)/%.tex
	@echo "Generating $@..."
	pdflatex -output-directory=$(DOC_DIR) $<

# Clean up generated files
clean:
	@echo "Cleaning up generated files..."
	rm -rf $(BIN_DIR)/*
	rm -f $(DOC_DIR)/*.pdf
	rm -f $(DOC_DIR)/*.aux $(DOC_DIR)/*.log
	rm -f $(PLOT_DIR)/plot
	rm -f $(PLOT_DIR)/data/*.txt $(PLOT_DIR)/figures/*.png
