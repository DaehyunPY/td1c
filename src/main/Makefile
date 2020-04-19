include ../../config/makeconfig

EXE1 = td1c.exe
#EXE1 = td1c_debug.exe

EXE2 = pes.exe
EXE3 = calprop.exe
EXE4 = test_field.exe
EXE5 = rtdisp.exe
EXE6 = read_print.exe

OBJ1 = \
td1c.o \
guess.o \
init.o \
tdse.o

OBJ2 = \
pes.o \
guess.o \

OBJ3 = \
calprop.o \
guess.o \
tdse.o \

OBJ4 = \
test_field.o

OBJ5 = \
rtdisp.o \
guess.o \
init.o \

OBJ6 = \
read_print.o \
guess.o \
tdse.o

#all: exe1 exe2 exe3 exe4 exe6
all: exe1

exe1: $(OBJ1)
	$(MPICXX) \
	../module/*.o \
	../io/*.o \
	../control/*.o \
	../mpi/*.o \
	../bas/*.o \
	../fedvr/*.o \
	../dvr/*.o \
	../sph/*.o \
	../tdcc/*.o \
	../ormas/*.o \
	../wfn/*.o \
	../hprod/*.o \
	../prop/*.o \
	../field/*.o \
	../gbas/*.o \
	../pes/*.o \
	../surff/*.o \
	../lapack/*.o \
	../util/*.o \
	../futil/*.o \
	$(OBJ1) \
	-I$(TD1C_INC) \
	$(LFLAG) \
	-o $(EXE1)
	@echo

exe2: $(OBJ2)
	$(MPICXX) \
	../module/*.o \
	../io/*.o \
	../control/*.o \
	../mpi/*.o \
	../bas/*.o \
	../fedvr/*.o \
	../dvr/*.o \
	../sph/*.o \
	../tdcc/*.o \
	../ormas/*.o \
	../wfn/*.o \
	../hprod/*.o \
	../prop/*.o \
	../field/*.o \
	../gbas/*.o \
	../pes/*.o \
	../surff/*.o \
	../lapack/*.o \
	../util/*.o \
	../futil/*.o \
	$(OBJ2) \
	-I$(TD1C_INC) \
	$(LFLAG) \
	-o $(EXE2)
	@echo

exe3: $(OBJ3)
	$(MPICXX) \
	../module/*.o \
	../io/*.o \
	../control/*.o \
	../mpi/*.o \
	../bas/*.o \
	../fedvr/*.o \
	../dvr/*.o \
	../sph/*.o \
	../tdcc/*.o \
	../ormas/*.o \
	../wfn/*.o \
	../hprod/*.o \
	../prop/*.o \
	../field/*.o \
	../gbas/*.o \
	../pes/*.o \
	../surff/*.o \
	../lapack/*.o \
	../util/*.o \
	../futil/*.o \
	$(OBJ3) \
	-I$(TD1C_INC) \
	$(LFLAG) \
	-o $(EXE3)
	@echo

exe4: $(OBJ4)
	$(MPICXX) \
	../module/*.o \
	../io/*.o \
	../control/*.o \
	../mpi/*.o \
	../field/*.o \
	$(OBJ4) \
	-I$(TD1C_INC) \
	$(LFLAG) \
	-o $(EXE4)
	@echo

exe5: $(OBJ5)
	$(MPICXX) \
	../module/*.o \
	../io/*.o \
	../control/*.o \
	../mpi/*.o \
	../bas/*.o \
	../fedvr/*.o \
	../dvr/*.o \
	../sph/*.o \
	../tdcc/*.o \
	../ormas/*.o \
	../wfn/*.o \
	../hprod/*.o \
	../prop/*.o \
	../field/*.o \
	../gbas/*.o \
	../pes/*.o \
	../lapack/*.o \
	../util/*.o \
	../futil/*.o \
	$(OBJ5) \
	-I$(TD1C_INC) \
	$(LFLAG) \
	-o $(EXE5)
	@echo

exe6: $(OBJ6)
	$(MPICXX) \
	../module/*.o \
	../io/*.o \
	../control/*.o \
	../mpi/*.o \
	../bas/*.o \
	../fedvr/*.o \
	../dvr/*.o \
	../sph/*.o \
	../tdcc/*.o \
	../ormas/*.o \
	../wfn/*.o \
	../hprod/*.o \
	../prop/*.o \
	../field/*.o \
	../gbas/*.o \
	../pes/*.o \
	../surff/*.o \
	../lapack/*.o \
	../util/*.o \
	../futil/*.o \
	$(OBJ6) \
	-I$(TD1C_INC) \
	$(LFLAG) \
	-o $(EXE6)
	@echo

include $(TD1C_CON)/make.objcomp
