# 指定编译器
CXX = g++
# 编译选项
CXXFLAGS = -g

# 自动匹配当前目录下所有以 Problem 开头的 .cpp 文件
PROBLEM_SRCS = $(wildcard Problem*.cpp)
PROBLEM_TARGETS = $(PROBLEM_SRCS:.cpp=)

# 报告的源文件（假设你用 LaTeX 写报告）
TEXFILE = report
TEX = xelatex
BIB = bibtex

# 默认规则，生成所有的目标可执行文件
all: $(PROBLEM_TARGETS)

# 生成每个可执行文件的规则（%.cpp -> %）
%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

# 编译并运行所有 Problem 的目标
run: $(PROBLEM_TARGETS)
	@for prog in $(PROBLEM_TARGETS); do \
		echo "Running $$prog..."; \
		./$$prog; \
		echo ""; \
	done

# 编译 LaTeX 报告生成 PDF

all: pdf

pdf: $(TEXFILE).pdf

report: $(TEXFILE).tex
	# First xelatex pass
	$(TEX) $(TEXFILE).tex
	# Run bibtex
	$(BIB) $(TEXFILE)
	# Second and third xelatex pass
	$(TEX) $(TEXFILE).tex
	$(TEX) $(TEXFILE).tex


# 清理所有生成的可执行文件和其他中间文件
clean:
	rm -f $(PROBLEM_TARGETS) *.out *.o *.aux *.bbl *.blg *.log *.out *.toc *.pdf *.run.xml *-blx.bib *.synctex.gz