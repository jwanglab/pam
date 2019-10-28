CC     = gcc
LIBS   = -lz -lm

OBJECTS = pam

all: $(OBJECTS)

pam:
	$(CC) $(CFLAGS) -o pam main.c $(LIBS)

.PHONY: clean
clean:
	-rm $(OBJECTS)
