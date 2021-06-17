F90       = mpif90
FCFLAGS = -O2 -fp-model=precise

BIN =  mm.ex pi.ex master.ex
OBJ = mm.o pi.o master.o

all: $(BIN)

%.o: %.f90
	$(F90) $(FCFLAGS) -c $^

mm.ex: mm.o
	$(F90) $(FCFLAGS) -o $@ mm.o

pi.ex: pi.o
	$(F90) $(FCFLAGS) -o $@ pi.o

master.ex: master.o
	$(F90) $(FCFLAGS) -o $@ master.o

clean:
	$(RM) $(BIN) $(OBJ)
