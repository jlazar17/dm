/** 

 @mainpage C++ Interface to Tauola 
 @brief Description of Tauola Interface in C++

 @authors Nadia Davidson, Gizo Nanava, Tomasz Przedzinski, Elzbieta Richter-Was, Zbigniew Was



 @section download1 New release

 The source code and documentation for release 1.0.7
 - <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/TAUOLA.1.0.7/Tauola_interface_design.1.0.7.pdf">Tauola_interface_design.1.0.7.pdf</a> full software documentation includes updates with respect to
the preprint version.
 - <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/TAUOLA.1.0.7/TAUOLA.1.0.7-LHC.tar.gz">TAUOLA source code for the LHC</a>
(or <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/TAUOLA.1.0.7/TAUOLA.1.0.7.tar.gz">TAUOLA source code with full TAUOLA FORTRAN</a>)
tarball and its <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/TAUOLA.1.0.7/svn_info_tauola.1.0.7.txt">revision info</a> SVN tag, tarball
creation date/time, etc. For updates  with respect to release 1.0 see <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/TAUOLA.1.0.7/changelog.1.0.7.txt">changelog.txt</a>.

- Note that LCG/Genser
<a href="http://sftweb.cern.ch/generators/">Generator
Services Subproject </a> distributes compiled, platform adapted
tar balls of our programs, in particular TAUOLA.1.0.2, TAUOLA.1.0.4 and TAUOLA.1.0.7.

 @section developement Developement version

 The source code and documentation are updated daily from the repository as well. The following files are provided for download of the developement version:
 - <a href="../Tauola_interface_design.pdf">Tauola_interface_design.pdf</a> full software documentation.
 - <a href="../TAUOLA.daily_temp.tar.gz">TAUOLA source code with full TAUOLA FORTRAN</a> tarball and its <a href="../svn_info_tauola.txt">revision info</a> SVN tag, tarball creation date/time, etc.
For updates  with respect to release 1.0 see <a href="changelog.txt">changelog.txt</a>. To remove code spurious for LHC C++ use, 
execute the script 'onlyLHC.sh' in the directory 'tauola-fortran'. Changes will be irreversible.

 @section download Older releases

 The source code and documentation for release 1.0. The following files are provided for download:
 - <a href="http://arxiv.org/abs/1002.0543">arXiv:1002.0543</a> full software documentation.
 - <a href="../TAUOLA.1.0.tar.gz">TAUOLA source code </a> tarball.

 The source code and documentation for release 1.0.2
 - <a href="../Tauola_interface_design.1.0.2.pdf">Tauola_interface_design.1.0.2.pdf</a> full software documentation includes updates with respect to 
the preprint version.
 - <a href="http://annapurna.ifj.edu.pl/~wasm/TAUOLA.1.0.2-LHC.tar.gz">TAUOLA source code for LHC</a>
(or <a href="../TAUOLA.1.0.2.tar.gz">TAUOLA source code with full TAUOLA FORTRAN</a>)
tarball and its <a href="../svn_info_tauola.1.0.2.txt">revision info</a> SVN tag, tarball
creation date/time, etc. For updates  with respect to release 1.0 see <a href="../changelog.1.0.2.txt">changelog.txt</a>.

 The source code and documentation for release 1.0.4
 - <a href="../Tauola_interface_design.1.0.4.pdf">Tauola_interface_design.1.0.4.pdf</a> full software documentation includes updates with respect to
the preprint version.
 - <a href="../TAUOLA.1.0.4-LHC.tar.gz">TAUOLA source code for the LHC</a>
(or <a href="../TAUOLA.1.0.4.tar.gz">TAUOLA source code with full TAUOLA FORTRAN</a>)
tarball and its <a href="../svn_info_tauola.1.0.4.txt">revision info</a> SVN tag, tarball
creation date/time, etc. For updates  with respect to release 1.0 see <a href="../changelog.1.0.4.txt">changelog.txt</a>.

 The source code and documentation for release 1.0.5
 - <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/TAUOLA.1.0.5/Tauola_interface_design.1.0.5.pdf">Tauola_interface_design.1.0.5.pdf</a> full software documentation includes updates with respect to
the preprint version.
 - <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/TAUOLA.1.0.5/TAUOLA.1.0.5-LHC.tar.gz">TAUOLA source code for the LHC</a>
(or <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/TAUOLA.1.0.5/TAUOLA.1.0.5.tar.gz">TAUOLA source code with full TAUOLA FORTRAN</a>)
tarball and its <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/TAUOLA.1.0.5/svn_info_tauola.1.0.5.txt">revision info</a> SVN tag, tarball
creation date/time, etc. For updates  with respect to release 1.0 see <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/TAUOLA.1.0.5/changelog.1.0.5.txt">changelog.txt</a>.

 @section intro Introduction/Status

 At present (since Feb 2 2010) the C++ interface functionality for TAUOLA is  complete.
 Longitudinal spin correlation are checked using <a href="http://mc-tester.web.cern.ch/MC-TESTER/">
 MC-TESTER</a> now. Transverse spin correlations and full functionality of TAUOLA  Universal Interface 
 are coded. Genuine electroweak corrections for processes mediated by Z/gamma are implemented as well.

 The tar file contains the C++ interface along with the source code for tauola itself (as available from <a href="http://wasm.web.cern.ch/wasm/goodies.html">old web page</a>; version Oct 11 2005).
 The developement version contains the latest source code for the interface from our subversion repository. Note that 
 revision numbers, dates and other info for this development version can be found in the revision information file listed above. 

 At present, the tar ball includes everything that is needed for installation. The user is advised to use  HepMC 2.04+ libraries as described below.

 @section setup Requirements

 For compilation, and to run the simple example, the interface requires:
 - <a href="http://lcgapp.cern.ch/project/simu/HepMC/">HepMC v2.04</a> or later.

 For further examples, one need to also install:
 - <a href="http://root.cern.ch/drupal/">ROOT v5.18</a> or later
 - <a href="http://home.thep.lu.se/~torbjorn/Pythia.html">PYTHIA 8.1</a> or later. PYTHIA must be compiled with HepMC 2 so that the PYTHIA library hepmcinterface exists.
 - <a href="http://mc-tester.web.cern.ch/MC-TESTER/">MC-TESTER v1.24</a> or later. Do not forget to compile the additional HepMC library libHepMCEvent as well.

 Examples of TAUOLA C++ combined with PHOTOS and its C++ interface are available from:
- <a href="http://www.ph.unimelb.edu.au/~ndavidson/photos/doxygen/index.html">
PHOTOS C++</a> official webpage.



 @section compile Configuration and Compilation

 In order to compile the TAUOLA C++ interface:
 - Execute './configure' with additional command line options:
   - '--with-hepmc=\<path\> ' provides the path to HepMC installation directory. One can set HEPMCLOCATION variable instead of using this directive. This options is required for interface to compile. To compile without HepMC use '--without-hepmc'
   - '--prefix=\<path\>' provides the installation path. The 'include' and 'lib' directories will be copied there if 'make install' is executed later. If none has been provided, the default directory for installation is '/usr/local'.
 - Execute 'make'
 - Optionally, execute 'make install' to copy files to the directory provided during configuration.

 After compiling the 'tauola-fortran' part, the TAUOLA C++ interface will be compiled and the '/lib' and '/include' directories will contain the appropriate library and include files.

 In order to compile the examples, enter the 'examples' directory, and:
 - execute './configure' to determine which examples can be compiled. Additional paths can be provided as command line options:
   - '--with-pythia8=\<path\>' provides the path to the Pythia8 installation directory. One can set the PYTHIALOCATION variable instead of using this directive. This path is required for all additional examples and tests.
   - '--with-mc-tester=\<path\>' provides the path to the MC-Tester installation directory (libHepMCEvent must be compiled as well, check MC-Tester documentation for more details). One can set the MCTESTERLOCATION variable instead of using this directive. This path is required for all additional examples and tests. It is assumed that using this option also implies that ROOT has already been installed (since it's required by MC-TESTER). The location of its binaries should be listed in the PATH variable.
 - execute 'make'

 Note that for examples working with PYTHIA 8.1, the PYTHIA8DATA global variable must be set (refer to instructions provided during configuration).
 Similarly, for examples in the examples/testing directory to work, the MCTESTERLOCATION global variable must be set.
 If neither PYTHIA nor MC-TESTER are present, only the simple example will be provided. The '/examples' directory will contain the compiled example files.

 If you prefer not to use configuration scripts,
 you may find <a href="README-NO-CONFIG.txt">README-NO-CONFIG.txt</a> useful.

 @section testing Testing

 In order to run some more specific tests both PYTHIA and MC-TESTER must be installed.
 - Compile TAUOLA C++ interface as well as examples.
 -  Check that the appropriate system variables are set: normally set by the script
 configure.paths.sh [.csh] (the configuration step mentions this script).
 -  Enter the /examples/testing directory. Modify test.inc if needed.
 - Enter the selected directory and execute 'make'.

 The appropriate .root files as well as .pdf files generated by MC-TESTER will be created inside the chosen directory. You can execute 'make clobber' to clean the directory. You can also execute 'make' inside the 'TAUOLA/examples/testing' directory to run all available tests one after another.

 @section sanc SANC
 Electroweak corrections may affect spin correlations between tau+ and tau- in a significant way. This is the case at high energies, 
far above WW threshold. At present we use the alpha scheme for electroweak calculations, this choice may not be optimal for the
cross section weight aiming at implementation of electroweak corrections. For this, the end scheme should be consistent with the choice 
of the host generator. On the other hand this is not the problem for spin correlations, our prime interest. In addition, if tables are used, transverse spin correlations for processes  mediated by Z and gamma are taken into account.

Changes to the tables generated by the SANC module can be implemented by modifying the interface located in the '/SANC' directory. Details regarding the structure and initialization variables of the interface are described in the documentation. In order to generate new tables:
 - execute 'make' in the '/SANC' directory to compile the library and tools.
 - Adjust the initialization in the SANC tables calculation if needed.
 - execute 'make tables' in the '/SANC' directory to generate the tables.
 - move tables to the directory from which your main program is executed.

 <hr>
 @section description Description of the code

 @image html tauola_interface_design.png "Design of Interface. Components are described in more detail below."
 
 @subsection outline Algorithm Outline

 The simplified version of the main program structure including the tauola decay routines would consist of:
 - initialize the Monte Carlo generator and all other generators / analysis tools.
 - set the appropriate parameters for the Tauola interface.
 - invoke the Tauola::initialize() routine.
 - for each event to be generated:
   - generate an event using any Monte Carlo generator
   - convert a generated event to the HepMC event record 
   - create a TauolaHepMCEvent from a HepMC event object
   - invoke decayTaus() for the TauolaHepMCEvent object. If a tau has already been decayed by a MC generator, use undecayTaus() first.
   - the HepMC event is now prepared and can be used in the analysis
 - finalize the analysis
 - print any summary output if needed and exit

 @subsection flow Flow of Control

 The following is a more detailed example point by point of how a tau, or set of taus, are decayed using the interface with the HepMC event record implementation.

 Starting from the main program:
 - At the initialization step, appropriate configuration routines via the 'Tauola' class methods can be executed.
 - The Tauola::initialize() method must be invoked, which calls the Initialize() methods in f_Init, which in turn calls the initilization routines of TAUOLA.

 In the main event loop, after the HepMC event is filled by a Monte Carlo generator:
 - An appropriate implementation of the abstract class TauolaEvent - TauolaHepMCEvent object must be created, taking a HepMC::GenEvent object as its parameter.
 - The method decayTaus() for the created TauolaHepMCEvent object must be invoked. If the Monte Carlo generator has already performed the decays of a tau, the undecayTaus() method can be executed just before the decayTaus() method.

 After executing decayTaus(), the following procedure takes place:
 - The HepMC event record is traversed and a list of stable taus in the event is created.
 - the taus from the list are paired with other decay products (e.g. tau neutrino or second tau). If the pair consists of two taus, the second tau is removed from the previous list.
 - TauolaParticlePairs are then created from these pairings.
 - For each pair, the spin density matrix is calculated using information about the production process.
 - The decayTauPairs() method is called for each pair. Within this method:
    - For each tau in the pair, decay() is called.
       - the decay method in TauolaParticle is responsible for adding itself
       (the tau) to the DecayList in the first position. It then calls
       decay() in f_Decay which in turn calls the dekay_() routine of TAUOLA.
       - dekay_() in TAUOLA generates a set of daughter particles with 
       momentum as though the tau was decayed in its rest frame.
    - The dekay_() also provides a polarimetric vector, which is used in combination with spin density matrix to calculate the spin weight. It is used to determine whether the decay is accepted or not.
    - If rejected, the list of decay products is cleared and the pair is decayed anew. This way unweighting of spin effects is performed.
    - Once accepted, for each tau in the pair the routine addDecayToEventRecord() is executed, which in turn calls TAUOLA dekay_() with an option set to write to the event record.
       - For each daughter TAUOLA calls filhep_() in the file f_FilHep.c.
       - filhep_() adds a new TauolaHepMCParticle to DecayList with all necessary information taken from TAUOLA.
       - The DecayList data structure translates the parent index given by TAUOLA to a pointer to the TauolaParticle(TauolaHepMCParticle) object.
       - The parents are set via the TauolaParticle methods (setMother, etc.) By setting the mothers it is automatically assigned a barcode and added to the HepMC event.
       - The particle is boosted by the tau's momentum.
       - The position of the vertex of the particle is modified accordingly to the tau lifetime.
 - Control is then returned to TauolaHepMCEvent which proceeds onto the next pair.
 - When all pairs are decayed, the TauolaHepMCEvent converts the HepMC::GenEvent to the units selected by the user (GeV/MeV mm/cm) and finalizes execution. The loop in the main program proceeds to the next event.

 @section event_rec Notes About the Event Record
 This section in future will  describe examples of events, particularly unusual ones,
 which the TAUOLA Interface  will nonetheless need to be able to handle, but which were not discussed at the time
of writing the documentation. At present this point is basically empty, except one example. Such type of events is already discussed in 
the documentation.
 - <a href="../pythia_event.pdf">example of a pythia event</a> with 
 non-conservation of momentum and tau->tau type vertices.

 @section to_do Things to do

See software documentation in pdf form. 


 <hr>
Last update March 20 2010. 
*/
