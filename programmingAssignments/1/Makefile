# 指定编译器
CXX = g++
# 编译选项
CXXFLAGS = -g

# 目标可执行文件
SRCS = B.cpp C.cpp D.cpp
TARGETS = $(SRCS:.cpp=)

# 默认规则，生成所有的目标可执行文件
all: $(TARGETS)

# 生成每个可执行文件的规则（%.cpp -> %）
%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

# 清理所有生成的可执行文件和其他中间文件
clean:
	rm -f $(TARGETS) *.o