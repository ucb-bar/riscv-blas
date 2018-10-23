include ../../make.inc

all: single double complex complex16
single:    xblat1s xblat2s xblat3s
double:    xblat1d xblat2d xblat3d
complex:   xblat1c xblat2c xblat3c
complex16: xblat1z xblat2z xblat3z

xblat1s: sblat1.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)
xblat1d: dblat1.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)
xblat1c: cblat1.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)
xblat1z: zblat1.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)

xblat2s: sblat2.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)
xblat2d: dblat2.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)
xblat2c: cblat2.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)
xblat2z: zblat2.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)

xblat3s: sblat3.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)
xblat3d: dblat3.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)
xblat3c: cblat3.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)
xblat3z: zblat3.o $(BLASLIB)
	$(LOADER) $(LOADOPTS) -o $@ $^ $(LOADLIBS)

run: all
	spike $(SPIKEARGS) pk xblat1s > sblat1.out
	spike $(SPIKEARGS) pk xblat1d > dblat1.out
	spike $(SPIKEARGS) pk xblat1c > cblat1.out
	spike $(SPIKEARGS) pk xblat1z > zblat1.out
	spike $(SPIKEARGS) pk xblat2s < sblat2.in
	spike $(SPIKEARGS) pk xblat2d < dblat2.in
	spike $(SPIKEARGS) pk xblat2c < cblat2.in
	spike $(SPIKEARGS) pk xblat2z < zblat2.in
	spike $(SPIKEARGS) pk xblat3s < sblat3.in
	spike $(SPIKEARGS) pk xblat3d < dblat3.in
	spike $(SPIKEARGS) pk xblat3c < cblat3.in
	spike $(SPIKEARGS) pk xblat3z < zblat3.in
	cat *.out | grep FAIL ; test $$? -eq 1

clean: cleanobj cleanexe cleantest
cleanobj:
	rm -f *.o
cleanexe:
	rm -f xblat*
cleantest:
	rm -f *.out core

.f.o:
	$(FORTRAN) $(OPTS) -c -o $@ $<
