# ====================================== #
# ============== Makefile ============== #
# ====================================== #

# ============ DIRECTORIES ============= #
SRC_DIR = src
OBJ_DIR = obj
EXE_DIR = exe
INC_DIR = include
DOC_DIR = doc/html

# ============= SOURCES FILES ========== #
FILES = $(notdir $(wildcard $(SRC_DIR)/*.cpp))

# ============= OUTPUT NAMES =========== #
PRODUCT   = ffast
ARCHIVE   = ffast

# ========== COMPILER & FLAGS ========== #

# Find if it is a Mac OS or Linux
OS := $(shell uname -s)

ifeq ($(OS),Linux)
    CC = g++
endif
ifeq ($(OS),Darwin)
    CC = clang++
endif

CPPFLAGS="-I/opt/local/include" 
LDFLAGS="-L/opt/local/lib"

CFLAGS = -Wall -O3 -std=c++11 -ffast-math

ifeq ($(OS),Darwin)
    CFLAGS += -stdlib=libc++
endif

OBJFLAGS = -lm -lfftw3

ifeq ($(OS),Linux)
    OBJFLAGS += -lrt
endif

CINCLUDE = -I $(INC_DIR)

# ======= SOURCES AND OBJECTS LIST ===== #
SRC = $(addprefix ${SRC_DIR}/, $(FILES))
OBJ = $(addprefix ${OBJ_DIR}/, $(addsuffix .o, $(basename $(FILES))))

# ========= DIRECTORY CREATION ========= #
MKDIR_P = mkdir -p

# ========== FORCING DELETION ========== #
RMFOR_P = rm -f
RMFOR_P_RECURS = $(RMFOR_P) -r

# ============= COMPILATION ============ #
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	@$(MKDIR_P) $(OBJ_DIR)
	$(CC) $(CFLAGS) $(CINCLUDE) $(CPPFLAGS) -c $< -o $@
$(EXE_DIR)/$(PRODUCT): $(OBJ)
	@$(MKDIR_P) $(EXE_DIR)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(OBJFLAGS)
	@echo ""
	@echo "The compilation is done."

# =========== OTHER COMMANDS =========== #
.PHONY:
	clean mrproper run tar doc help

clean:
	@$(RMFOR_P) $(OBJ)
	@$(RMFOR_P) ${EXE_DIR}/${PRODUCT}
	@$(RMFOR_P_RECURS) -r ${DOC_DIR}
	@echo "The compilation and doc files have been deleted."

mrproper: clean
	@$(RMFOR_P) ${EXE_DIR}/${PRODUCT}
	@echo "The final executable has been deleted."

tar:
	@tar -cvf $(addsuffix .tar, ${ARCHIVE}) Makefile readme.md Paper.pdf ${SRC_DIR} ${INC_DIR} > /dev/null
	@gzip -f9 $(addsuffix .tar, ${ARCHIVE})
	@echo "The archive has been generated."

doc: Doxyfile $(SRC)
	@doxygen

help:
	@echo ""
	@echo "List of available commands :"
	@echo "    - clean : Delete all the file created during the compilation."
	@echo "    - mrproper : Delete only the executable file"
	@echo "    - tar : Create an archive with all the files needed for compilation."
	@echo ""
