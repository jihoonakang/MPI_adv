F90       = mpif90
FCFLAGS = -O2
LFLAGS = -fPIC

BIN =  mm.ex pi.ex master.ex
OBJ = mm.o pi.o master.o

all: $(BIN)

%.o: %.f90
	$(F90) $(FCFLAGS) -c $^

mm.ex: mm.o
	$(F90) $(FCFLAGS) $(LFLAGS) -o $@ mm.o

pi.ex: pi.o
	$(F90) $(FCFLAGS) $(LFLAGS) -o $@ pi.o

master.ex: master.o
	$(F90) $(FCFLAGS) $(LFLAGS) -o $@ master.o

clean:
	$(RM) $(BIN) $(OBJ)
