PAGES_DIR = doxygenPages

default: Doxyfile mainpage.h 
	doxygen Doxyfile
	@echo "******************************************************"
	@echo " Please load the doxygen documentation file           "
	@echo " documentation/html/index.html in your web browser    "	
	@echo "******************************************************"

clean:
	rm -rf html/
	rm -f *~
