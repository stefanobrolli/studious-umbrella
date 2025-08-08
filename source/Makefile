# Compiler and flags
CC = g++
CFLAGS = -O3 -ffast-math -march=native -flto -MMD -MP

# List of source files 
SRC = Main.cpp Basis.cpp Props.cpp SelfEn.cpp Functions.cpp

# Object files (each .cpp file is compiled to a corresponding .o file)
OBJ = $(SRC:.cpp=.o)

# Dependency files (.d files are generated automatically)
DEPS = $(OBJ:.o=.d)

# Target for the final executable
TARGET = RichModSCGF

# Rule to build the executable from object files
$(TARGET): $(OBJ)
	$(CC) -o $@ $^ 

# Rule for compiling each .cpp file into an .o file, and generate dependencies
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Include dependency files if they exist
-include $(DEPS)

# Phony target to clean up generated files
.PHONY: clean
clean:
	rm -f $(OBJ) $(DEPS) $(TARGET)



