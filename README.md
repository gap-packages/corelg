# CoReLG

CoReLG is a GAP4 package. Its main objective is to provide
functionality for computing with real (semi-)simple Lie algebras.


## Installation

To install the package CoReLG move the file `corelg-XX.tar.gz`
(or any other archive containing it) into the `pkg` directory.
Usually this will be the `pkg` subdirectory in your GAP4 installation.
However, it is also possible to have a `pkg` subdirectory in a 
different place, see the section `Installing GAP Packages` of the 
GAP4 reference manual for more information.
Then simply unpack `corelg-XX.tar.gz` and your installation is
complete.
In GAP issue 

    gap> LoadPackage( "corelg" );

             
## Documentation

The manual of CoReLG is contained in the `doc` directory. The different 
versions of the manual (pdf, dvi, html) can be compiled by doing

    gap> corelg.makeManual( PATHDOC, PATHGAP );

where PATHDOC is the path to the /pkg/corelg/doc directory, and
PATHGAP is the path to the gap root directory.
