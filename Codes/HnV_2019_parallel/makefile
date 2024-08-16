CC = gcc
CFLAGS = -Wall -O2 -fopenmp -lm -g
LDFLAGS = -lm -shared
SRC_DIR = .
BUILD_DIR = build
BIN_DIR = .

TARGET = LBM.so
SOURCES = $(wildcard $(SRC_DIR)/*.c)
OBJECTS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SOURCES))

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(BIN_DIR)/$@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) -fPIC -c -o $@ $<

clean:
	rm -rf $(BUILD_DIR)/*.o $(BIN_DIR)/$(TARGET)

.PHONY: clean

