CC = g++
CFLAGS = -W -Wall -O2 -std=c++14 -I./include
CL_FLAGS =
TARGET = test.a

UNAME_S = $(shell uname -s)

ifeq ($(UNAME_S), Linux)
	CL_FLAGS = -L/usr/local/cuda/lib64 -lOpenCL
endif
ifeq ($(UNAME_S), Darwin)
    CL_FLAGS = -framework OpenCL
endif

DIRS = objs/ objs/clmpz/ objs/clmpr objs/ckks/

CLMP_OBJS = objs/clmpr/rand_bits_ui.o objs/clmpr/rand_bnd_ui.o \
	objs/clmpz/set.o objs/clmpz/set_ui.o objs/clmpz/set_ul.o objs/clmpz/set_i.o objs/clmpz/set_li.o \
	objs/clmpz/add.o objs/clmpz/add_ui.o objs/clmpz/add_ul.o \
	objs/clmpz/bit.o objs/clmpz/bits.o \
	objs/clmpz/shift_left_ui.o objs/clmpz/shift_right_ui.o \
	objs/clmpz/sign.o objs/clmpz/cmp.o objs/clmpz/is_zero.o \
	objs/clmpz/bit_and.o objs/clmpz/bit_or.o objs/clmpz/bit_xor.o \
	objs/clmpz/neg.o \
	objs/clmpz/sub.o objs/clmpz/sub_ui.o objs/clmpz/sub_ul.o \
	objs/clmpz/mul.o objs/clmpz/mul_ui.o objs/clmpz/mul_ul.o\
	objs/clmpz/div.o objs/clmpz/div_ui.o objs/clmpz/div_ul.o \
	objs/clmpz/mod.o objs/clmpz/mod_ui.o objs/clmpz/mod_ul.o \
	objs/clmpz/rand_bits.o \
	objs/clmpz/to_string.o objs/clmpz/to_ui.o objs/clmpz/to_ul.o

CKKS_OBJS = objs/ckks/context.o objs/ckks/scheme.o \
	objs/ckks/key.o objs/ckks/secretkey.o \
	objs/ckks/plaintext.o objs/ckks/ciphertext.o


TARGET = test.a
TARGET_SRC = run/test.cpp

MKDIR_P = mkdir -p

UTIL_DIR = .../src/utility
UTIL_OBJ_DIR = .../objs/utility
UTIL_OBJ = $(UTIL_DIR)/%.o

all: $(DIRS) $(CLMP_OBJS) $(CKKS_OBJS) $(TARGET)

$(DIRS):
	$(MKDIR_P) $@

$(TARGET): $(TARGET_SRC) $(CLMP_OBJS) $(CKKS_OBJS)
	$(CC)  $(CFLAGS)  -o $@ $^ $(CL_FLAGS)

objs/clmpr/rand_bits_ui.o: src/clmpr/rand_bits_ui.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpr/rand_bnd_ui.o: src/clmpr/rand_bnd_ui.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/set.o: src/clmpz/set.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/set_ui.o: src/clmpz/set_ui.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/set_ul.o: src/clmpz/set_ul.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/set_i.o: src/clmpz/set_i.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/set_li.o: src/clmpz/set_li.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/add.o: src/clmpz/add.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/add_ui.o: src/clmpz/add_ui.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/add_ul.o: src/clmpz/add_ul.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/bit.o: src/clmpz/bit.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/bits.o: src/clmpz/bits.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/bit_and.o: src/clmpz/bit_and.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/bit_or.o: src/clmpz/bit_or.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/bit_xor.o: src/clmpz/bit_xor.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/shift_left_ui.o: src/clmpz/shift_left_ui.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/shift_right_ui.o: src/clmpz/shift_right_ui.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/sign.o: src/clmpz/sign.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/cmp.o: src/clmpz/cmp.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/is_zero.o: src/clmpz/is_zero.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/neg.o: src/clmpz/neg.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/sub.o: src/clmpz/sub.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/sub_ui.o: src/clmpz/sub_ui.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/sub_ul.o: src/clmpz/sub_ul.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/mul.o: src/clmpz/mul.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/mul_ui.o: src/clmpz/mul_ui.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/mul_ul.o: src/clmpz/mul_ul.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/div.o: src/clmpz/div.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/div_ui.o: src/clmpz/div_ui.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/div_ul.o: src/clmpz/div_ul.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/mod.o: src/clmpz/mod.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/mod_ui.o: src/clmpz/mod_ui.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/mod_ul.o: src/clmpz/mod_ul.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/rand_bits.o: src/clmpz/rand_bits.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/to_string.o: src/clmpz/to_string.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/to_ui.o: src/clmpz/to_ui.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/clmpz/to_ul.o: src/clmpz/to_ul.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/ckks/context.o: src/ckks/context.cpp
	$(CC) $(CFLAGS)  -c -o $@ $^

objs/ckks/key.o: src/ckks/key.cpp
	$(CC) $(CFLAGS)  -c -o $@ $^

objs/ckks/secretkey.o: src/ckks/secretkey.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/ckks/scheme.o: src/ckks/scheme.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/ckks/plaintext.o: src/ckks/plaintext.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

objs/ckks/ciphertext.o: src/ckks/ciphertext.cpp
	$(CC) $(CFLAGS) -c -o $@ $^

clean:
	rm -rf $(OBJS) $(TARGET) $(DIRS)
