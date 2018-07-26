############## Make tests ###########################
include ../test.inc

CONF=pythia_H.conf

testAllH: tester.H+.pdf tester.H-.pdf

## W -> tau nu, with tau -> pi nu, MC-TESTER using mass^2
tester.H+.pdf: mc-tester.H+.root
	$(COMPARE) mc-tester.H+.original.root mc-tester.H+.root
	mv tester.pdf tester.H+.pdf
	$(VIEWER) tester.H+.pdf &

mc-tester.H+.root: $(PROGRAM) $(CONF)
	cp SETUP_H+.C SETUP.C
	$(PROGRAM) $(CONF) 1 $(EVENTS) $(PI_MODE)
	mv mc-tester.root mc-tester.H+.root
	rm -rf SETUP.C

## W -> tau nu, with tau -> pi nu, MC-TESTER using mass^2
tester.H-.pdf: mc-tester.H-.root
	$(COMPARE) mc-tester.H-.original.root mc-tester.H-.root
	mv tester.pdf tester.H-.pdf
	$(VIEWER) tester.H-.pdf &

mc-tester.H-.root: $(PROGRAM) $(CONF)
	cp SETUP_H-.C SETUP.C
	$(PROGRAM) $(CONF) 1 $(EVENTS) $(PI_MODE)
	mv mc-tester.root mc-tester.H-.root
	rm -rf SETUP.C

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.H+.root 
	rm -f mc-tester.H-.root

