############## Make tests ###########################
include ../test.inc

CONF=pythia_W.conf

testAllW: tester.W+.pdf tester.W-.pdf

## W -> tau nu, with tau -> pi nu, MC-TESTER using mass^2
tester.W+.pdf: mc-tester.W+.root
	$(COMPARE) mc-tester.W+.original.root mc-tester.W+.root
	mv tester.pdf tester.W+.pdf
	$(VIEWER) tester.W+.pdf &

mc-tester.W+.root: $(PROGRAM) $(CONF)
	cp SETUP_W+.C SETUP.C
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(PI_MODE)
	mv mc-tester.root mc-tester.W+.root
	rm -rf SETUP.C

## W -> tau nu, with tau -> pi nu, MC-TESTER using mass^2
 tester.W-.pdf: mc-tester.W-.root
	$(COMPARE) mc-tester.W-.original.root mc-tester.W-.root
	mv tester.pdf tester.W-.pdf
	$(VIEWER) tester.W-.pdf &

mc-tester.W-.root: $(PROGRAM) $(CONF)
	cp SETUP_W-.C SETUP.C
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(PI_MODE)
	mv mc-tester.root mc-tester.W-.root
	rm -rf SETUP.C

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.W+.root 
	rm -f mc-tester.W-.root

