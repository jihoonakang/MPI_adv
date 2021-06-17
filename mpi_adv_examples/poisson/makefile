.SUFFIXES : .c .o 
CC = gcc
OBJS = poisson.o 

CFLAGS= -O2
LIBS	:= -lm

TARGET =  poisson.e

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf $(TARGET) $(OBJS) *.o

new:
	touch $(SRCS)

