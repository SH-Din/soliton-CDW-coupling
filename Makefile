zCXX = D:/ProgramFiles/msys64/mingw64/bin/g++.exe
CXXFLAGS = -std=c++17 -O3 -fopenmp \
	-Iinclude \
	-ID:/ScienceResearch/Eigen/eigen-3.3.9 \
	-ID:/ScienceResearch/json_reader/json/include

SRCDIR = source
OBJDIR = out

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))
TARGET  = $(OBJDIR)\free_energy_scanning.exe

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	powershell -Command "Remove-Item -Path '$(OBJDIR)\*.o','$(TARGET)' -Force -ErrorAction SilentlyContinue"