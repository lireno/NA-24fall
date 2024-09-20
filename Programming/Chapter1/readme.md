# Numerical Analysis Programming Homework

## Overview
This repository contains the implementation of numerical methods (Bisection, Newton's, and Secant methods) for solving various mathematical problems as part of the Numerical Analysis programming homework.

### Files:
- **EquationSolver.hpp**: Defines the base class for all root-finding algorithms.
- **Function.hpp**: Provides function wrappers for easy evaluation of mathematical expressions.
- **ProblemB.cpp**: Implements and tests the bisection method for specific functions.
- **ProblemC.cpp**: Implements and tests Newton's method for specific functions.
- **ProblemD.cpp**: Implements and tests the secant method for specific functions.
- **ProblemE.cpp**: Solves the trough problem using the Bisection, Newton's, and Secant methods.
- **ProblemF.cpp**: Solves the nose-in failure angle problem using Newton's and Secant methods.
- **Report.tex**: LaTeX source file for generating the report.
- **Makefile**: Automates the build and execution process.

## How to Use the Makefile

### 1. Compile and Run All Problem Code

To compile and run all the problem files, simply use the following command:

```bash
make run
```

This command compiles all the problem source files (ProblemB.cpp, ProblemC.cpp, ProblemD.cpp, ProblemE.cpp, ProblemF.cpp) and runs the resulting executables. The output for each problem will be displayed in the terminal.

### 2. Compile the Report and Generate a PDF

To compile the LaTeX report and generate a PDF file, run:

```bash
make report
```

This command compiles the LaTeX file `Report.tex` and generates `Report.pdf`, which contains all the necessary documentation for the assignment.

### 3. Clean the Directory

To clean up all object files, executables, and temporary files generated during compilation, use:

```bash
make clean
```

This will remove all the compiled files, executables, and temporary LaTeX files, leaving the directory clean.

### 4. Modify `ProblemB.cpp` and Handle Errors

In the source file `ProblemB.cpp`, there's a section of code that is currently commented out:

```cpp
// solve_f4();
```

To observe the error caused by this code, navigate to line 68 in `ProblemB.cpp`, remove the comment `//` before `solve_f4();`, and save the file. Then, recompile `ProblemB.cpp` and run the executable using the following commands:

```bash
make ProblemB
```

The program will attempt to solve a function that may have discontinuities or other issues, causing an error. The program should output an error message indicating the problem (e.g., function discontinuity).
