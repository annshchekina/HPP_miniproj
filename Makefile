CC = gcc
LD = gcc
CFLAGS = -Wall -std=c99
LDFLAGS = 
LIBS = -lm
RM = /bin/rm -f

OBJS = decrypt.o file_utils.o utils.o funcs.o timings.o

EXECUTABLE = decrypt

all:$(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXECUTABLE)

decrypt.o: decrypt.c file_utils.h utils.h funcs.h timings.h
	$(CC) $(CFLAGS) -c decrypt.c

file_utils.o: file_utils.c file_utils.h
	$(CC) $(CFLAGS) -c file_utils.c

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) -c utils.c

funcs.o: funcs.c funcs.h
	$(CC) $(CFLAGS) -c funcs.c

timings.o: timings.c timings.h
	$(CC) $(CFLAGS) -c timings.c

clean:
	$(RM) $(EXECUTABLE) $(OBJS) *~

