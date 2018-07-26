############## Make tests ###########################
include ../test.inc

CONF=pythia_H.conf

testAllH: tester.H.Pi.mass1.pdf  tester.H.Pi.mass2.pdf \
	tester.H.Rho.pdf tester.HA.Pi.pdf tester.HA.Rho.pdf

## Higgs -> tau tau, with tau -> pi nu, MC-TESTER using mass^1
tester.H.Pi.mass1.pdf: mc-tester.H.Pi.mass1.root
	$(COMPARE) mc-tester.H.Pi.mass1.original.root mc-tester.H.Pi.mass1.root
	mv tester.pdf tester.H.Pi.mass1.pdf
	$(VIEWER) tester.H.Pi.mass1.pdf &

mc-tester.H.Pi.mass1.root: $(PROGRAM) $(CONF)
	cp SETUP_H_PI_1.C SETUP.C
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(PI_MODE)
	mv mc-tester.root mc-tester.H.Pi.mass1.root
	rm -rf SETUP.C

## Higgs -> tau tau, with tau -> pi nu, MC-TESTER using mass^2
tester.H.Pi.mass2.pdf: mc-tester.H.Pi.mass2.root
	$(COMPARE) mc-tester.H.Pi.mass2.original.root mc-tester.H.Pi.mass2.root
	mv tester.pdf tester.H.Pi.mass2.pdf
	$(VIEWER) tester.H.Pi.mass2.pdf &

mc-tester.H.Pi.mass2.root: $(PROGRAM) $(CONF)
	cp SETUP_H_PI_2.C SETUP.C
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(PI_MODE)
	mv mc-tester.root mc-tester.H.Pi.mass2.root
	rm -rf SETUP.C

## Higgs -> tau tau, with tau -> rho nu
tester.H.Rho.pdf: mc-tester.H.Rho.root
	$(COMPARE) mc-tester.H.Rho.original.root mc-tester.H.Rho.root
	mv tester.pdf tester.H.Rho.pdf
	$(VIEWER) tester.H.Rho.pdf &

mc-tester.H.Rho.root: $(PROGRAM) $(CONF)
	cp SETUP_H_RHO.C SETUP.C
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(RHO_MODE)
	mv mc-tester.root mc-tester.H.Rho.root
	rm -rf SETUP.C

##################
tester.HA.Pi.pdf: mc-tester.HA.Pi.root
	$(COMPARE) mc-tester.HA.Pi.original.root mc-tester.HA.Pi.root
	mv tester.pdf tester.HA.Pi.pdf
	$(VIEWER) tester.HA.Pi.pdf &

mc-tester.HA.Pi.root: $(PROGRAM) $(CONF)
	cp SETUP_H_PI_1.C SETUP.C
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(PI_MODE) 0.7853
	mv mc-tester.root mc-tester.HA.Pi.root
	rm -rf SETUP.C

################
tester.HA.Rho.pdf: mc-tester.HA.Rho.root
	$(COMPARE) mc-tester.HA.Rho.original.root mc-tester.HA.Rho.root
	mv tester.pdf tester.HA.Rho.pdf
	$(VIEWER) tester.HA.Rho.pdf &

mc-tester.HA.Rho.root: $(PROGRAM) $(CONF)
	cp SETUP_H_RHO.C SETUP.C
	$(PROGRAM) $(CONF) 0 $(EVENTS) $(RHO_MODE) 0.7853
	mv mc-tester.root mc-tester.HA.Rho.root
	rm -rf SETUP.C

clean:
	rm -rf *~ booklet* tester* mc-results* 
	rm -rf MC-TESTER.DAT *_C.so *_C.d

clobber: clean
	rm -f mc-tester.H.Pi.mass1.root 
	rm -f mc-tester.H.Pi.mass2.root
	rm -f mc-tester.H.Rho.root 
	rm -f mc-tester.HA.Pi.root
	rm -f mc-tester.HA.Rho.root
