############## Make tests ###########################
include ../test.inc

CONF=pythia_Z.conf

testAllZ: tester.Z.Pi.mass1.pdf tester.Z.Pi.mass2.pdf
#	tester.Z.Rho.pdf

## Z -> tau tau, with tau -> pi nu, MC-TESTER using mass^1
tester.Z.Pi.mass1.pdf: mc-tester.Z.Pi.mass1.root
	$(COMPARE) mc-tester.Z.Pi.mass1.original.root mc-tester.Z.Pi.mass1.root
	mv tester.pdf tester.Z.Pi.mass1.pdf
	$(VIEWER) tester.Z.Pi.mass1.pdf &

mc-tester.Z.Pi.mass1.root: $(PROGRAM) $(CONF)
	cp SETUP_Z_PI_1.C SETUP.C
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(PI_MODE)
	mv mc-tester.root mc-tester.Z.Pi.mass1.root
	rm -rf SETUP.C

## Z -> tau tau, with tau -> pi nu, MC-TESTER using mass^2
tester.Z.Pi.mass2.pdf: mc-tester.Z.Pi.mass2.root
	$(COMPARE) mc-tester.Z.Pi.mass2.original.root mc-tester.Z.Pi.mass2.root
	mv tester.pdf tester.Z.Pi.mass2.pdf
	$(VIEWER) tester.Z.Pi.mass2.pdf &

mc-tester.Z.Pi.mass2.root: $(PROGRAM) $(CONF)
	cp SETUP_Z_PI_2.C SETUP.C
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(PI_MODE)
	mv mc-tester.root mc-tester.Z.Pi.mass2.root
	rm -rf SETUP.C

## Uncomment when we have Z transverse spin included
## Z -> tau tau, with tau -> rho nu 
#tester.Z.Rho.pdf: mc-tester.Z.Rho.root
#	$(COMPARE) mc-tester.Z.Rho.original.root mc-tester.Z.Rho.root
#	mv tester.pdf tester.Z.Rho.pdf
#	$(VIEWER) tester.Z.Rho.pdf &

#mc-tester.Z.Rho.root: $(PROGRAM) $(CONF)
#	cp SETUP_Z_RHO.C SETUP.C
#	$(PROGRAM) $(CONF) 0 $(EVENTS) $(RHO_MODE)
#	mv mc-tester.root mc-tester.Z.Rho.root
#	rm -rf SETUP.C
#	cp mc-tester.Z.Rho.root mc-tester.Z.Rho.original.root
#remove this last line once the benchmark file has been made

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.Z.Pi.mass1.root 
	rm -f mc-tester.Z.Pi.mass2.root
#	rm -f mc-tester.Z.Rho.root 
