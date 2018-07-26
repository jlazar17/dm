############## Make tests ###########################
include ../test.inc

CONF=pythia_Z.conf

testAllTaus: tester.tau+.pdf tester.tau-.pdf

## W -> tau nu, with tau -> pi nu, MC-TESTER using mass^2
tester.tau+.pdf: mc-tester.tau+.root
	$(COMPARE) mc-tester.tau+.original.root mc-tester.tau+.root
	mv tester.pdf tester.tau+.pdf
	$(VIEWER) tester.tau+.pdf &

mc-tester.tau+.root: $(PROGRAM) $(CONF)
	cp SETUP_tau+.C SETUP.C
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(ALL_MODE)
	mv mc-tester.root mc-tester.tau+.root
	rm -rf SETUP.C

## W -> tau nu, with tau -> pi nu, MC-TESTER using mass^2
tester.tau-.pdf: mc-tester.tau-.root
	$(COMPARE) mc-tester.tau-.original.root mc-tester.tau-.root
	mv tester.pdf tester.tau-.pdf
	$(VIEWER) tester.tau-.pdf &

mc-tester.tau-.root: $(PROGRAM) $(CONF)
	cp SETUP_tau-.C SETUP.C
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(ALL_MODE)
	mv mc-tester.root mc-tester.tau-.root
	rm -rf SETUP.C

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.tau+.root 
	rm -f mc-tester.tau-.root