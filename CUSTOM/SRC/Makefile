include ../../make.inc
all: custom-utils.o

clean: cleanobj cleanlib

cleanobj:
	rm -f *.o
cleanlib:

.f.o:
	$(FORTRAN) $(OPTS) -c -o $@ $<
.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<
