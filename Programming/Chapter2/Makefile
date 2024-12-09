# 使用 g++ 编译器
CXX = g++
# 编译选项
CXXFLAGS = -g -I./include  # 增加头文件路径

# 源文件目录和目标文件目录
SRCDIR = ./src
INCLUDEDIR = ./include
BINDIR = ./bin

# 搜索所有 src 文件夹下的 Problem*.cpp 文件
PROBLEM_SRCS = $(wildcard $(SRCDIR)/Problem*.cpp)
# 将 src 下的 cpp 文件名替换为 bin 下的可执行文件名
PROBLEM_TARGETS = $(patsubst $(SRCDIR)/Problem%.cpp, $(BINDIR)/Problem%, $(PROBLEM_SRCS))

# LaTeX 文件相关
TEXFILE = report
TEX = pdflatex
BIB = bibtex

all: $(PROBLEM_TARGETS) report

%: $(SRCDIR)/Problem%.cpp
	@mkdir -p $(BINDIR)  # 如果 bin 文件夹不存在，则创建
	$(CXX) $(CXXFLAGS) -o $(BINDIR)/Problem$@ $<
	$(BINDIR)/Problem$@

$(BINDIR)/Problem%: $(SRCDIR)/Problem%.cpp
	@mkdir -p $(BINDIR)  # 如果 bin 文件夹不存在，则创建
	$(CXX) $(CXXFLAGS) -o $@ $<
	./$@

# 运行所有 Problem 的目标
run: $(PROBLEM_TARGETS)
	@for prog in $(PROBLEM_TARGETS); do \
		echo "Running $$prog..."; \
		./$$prog; \
		echo ""; \
	done

report: $(TEXFILE).tex
	# First pdflatex pass
	$(TEX) $(TEXFILE).tex
	# Run bibtex
	$(BIB) $(TEXFILE)
	# Second and third pdflatex pass
	$(TEX) $(TEXFILE).tex
	$(TEX) $(TEXFILE).tex
	rm -f $(TEXFILE).aux $(TEXFILE).bbl $(TEXFILE).blg $(TEXFILE).log $(TEXFILE).out $(TEXFILE).toc $(TEXFILE).run.xml $(TEXFILE)-blx.bib $(TEXFILE).synctex.gz

clean:
	rm -rf $(BINDIR)/*.o *.aux *.bbl *.blg *.log *.out *.toc *.pdf *.run.xml *-blx.bib *.synctex.gz
