

SRC_PATH = ./src

BIN_DIR = ./bin
OBJ_DIR = ./obj

WFLAGS = -Wall -pedantic -Wshadow -Winline






ifeq ($(mode),debug)
   ADD_FLAGS = -g
else
   mode = release
   ADD_FLAGS = -O3
endif

ALL_FLAGS = $(WFLAGS) $(ADD_FLAGS)

ifeq ($(parallel),y)
	ALL_FLAGS = $(WFLAGS) $(ADD_FLAGS) -fopenmp
else
	ALL_FLAGS = $(WFLAGS) $(ADD_FLAGS)
endif


CC =g++
EXECUTEABLE = strike_pdbcontactsgenerator

all: $(BIN_DIR)/$(EXECUTEABLE)



CLASS_LIB = $(SRC_PATH)/classes
UTIL_LIB = $(SRC_PATH)/util
ALN_METHOD_LIB = $(SRC_PATH)/aln_methods


CLASS_OBJ = PDB.class_o Contacts.class_o Alignment.class_o Matrices.class_o
# RestServices.class_o
OBJ = main.o
UTIL_OBJ = aligning.util_o
# xml.util_o


# TEST_OBJ = test.o




_ALL_OBJ = $(CLASS_OBJ) $(UTIL_OBJ) $(OBJ)

ALL_OBJ = $(patsubst %,$(OBJ_DIR)/%,$(_ALL_OBJ))




$(BIN_DIR)/$(EXECUTEABLE): $(ALL_OBJ)
	$(CC) $(ALL_FLAGS) $(ALL_OBJ) -o $(BIN_DIR)/$(EXECUTEABLE)
# 	-lcurl
# 	cp $(BIN_DIR)/$(EXECUTEABLE) $(HOME)/bin


$(OBJ_DIR)/%.class_o: $(CLASS_LIB)/%.cpp $(CLASS_LIB)/%.h
	$(CC) $(ALL_FLAGS) -c $< -o $@

$(OBJ_DIR)/%.aln_o: $(ALN_METHOD_LIB)/%.cpp
	$(CC) $(ALL_FLAGS) -c $< -o $@

$(OBJ_DIR)/%.util_o: $(UTIL_LIB)/%.cpp
	$(CC) $(ALL_FLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_PATH)/%.cpp
	$(CC) $(ALL_FLAGS) -c $< -o $@


# pdb_test: $(ALL_OBJ)
# 	$(CC) $(ALL_FLAGS) -c $< -o $@
# 	$(CC) $(ALL_FLAGS) $(ALL_OBJ) -o $(BIN_DIR)/pdb_test -lcurl



.PHONY: doc readme

doc:
	doxygen $(DOC_DIR)/doc.doxygen

readme:
	doxygen $(DOC_DIR)/readme.doxygen

.PHONY: clean
clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(OBJ_DIR)/*.class_o
	rm -f $(OBJ_DIR)/*.aln_o
	rm -f $(OBJ_DIR)/*.util_o
	rm -f $(BIN_DIR)/$(EXECUTEABLE)
	rm -rf $(DOC_DIR)/html $(DOC_DIR)/latex





