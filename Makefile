CC = gcc
# CFLAGS = -Iinclude -Wall -Wextra -Werror -fPIC
CFLAGS = -Iinclude -Wall -Wextra -fPIC
# SRC_DIR = src
# SRC = $(wildcard $(SRC_DIR)/*.c)
SRC = src/mat_utils.c src/matrix.c
TEST_DIR = test
TEST_SRC = $(wildcard $(TEST_DIR)/*.c)
OBJ = $(SRC:.c=.o)
MAIN = main.c
MAIN_OBJ = $(MAIN:.c=.o)
TEST_OBJ = $(TEST_SRC:.c=.o)
TARGET = main
LIB = libmatrix.so
TEST_TARGETS = $(TEST_OBJ:.o=)

all: $(TARGET) $(TEST_TARGETS) $(LIB)

# Compile source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Build the shared library
$(LIB): $(OBJ)
	$(CC) -shared -o $(LIB) $(OBJ)

# Link the main program
$(TARGET): $(OBJ) $(MAIN_OBJ)
	$(CC) $(OBJ) $(MAIN_OBJ) -o $(TARGET)

# test program
test/%: test/%.o $(OBJ)
	$(CC) $^ -o $@

clean:
	rm -f $(OBJ) $(MAIN_OBJ) $(TEST_OBJ) $(TARGET) $(TEST_TARGETS) $(LIB)

.PHONY: all clean