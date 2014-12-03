#####################################################################################
#
#  PackageInfo.g                    Heiko Dietrich, Paolo Faccin, and Willem de Graaf 
#
#
# The package CoReLG is free software; you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version. 

# Details may have to be corrected...


SetPackageInfo( rec(
PackageName := "CoReLG",
Subtitle := "computation with real Lie groups",        
Version := "1.20",
Date := "02/12/2014",
ArchiveURL := Concatenation("http://users.monash.edu/~heikod/corelg/corelg-",~.Version),
ArchiveFormats := ".tar.gz",
Persons := [

  rec(
  LastName := "Dietrich",
  FirstNames := "Heiko",
  IsAuthor := true,
  IsMaintainer := true,
  Email := "heiko.dietrich@monash.edu",
  WWWHome := "http://users.monash.edu.au/~heikod/",
  Place := "Melbourne",
  Institution := "School of Mathematical Sciences, Monash University"
  ),

  rec(
  LastName := "Faccin",
  FirstNames := "Paolo",
  IsAuthor := true,
  IsMaintainer := true,
  Email := "faccin@science.unitn.it",
  Place := "Trento",
  Institution := "Dipartimento di Matematica, University of Trento"
  ),

  rec(
  LastName := "de Graaf",
  FirstNames := "Willem Adriaan",
  IsAuthor := true,
  IsMaintainer := true,
  Email := "degraaf@science.unitn.it",
  WWWHome := "http://www.science.unitn.it/~degraaf",
  Place := "Trento",
  Institution := "Dipartimento di Matematica, University of Trento"
  )

],
Status := "accepted",
CommunicatedBy := "Bettina Eick (Braunschweig)",
AcceptDate := "01/2014",
PackageDoc := rec( BookName  := "CoReLG" ,  
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Computing with real Lie groups",
  Autoload  := false
),
README_URL := 
  "http://users.monash.edu/~heikod/corelg/README",
PackageInfoURL := 
  "http://users.monash.edu/~heikod/corelg/PackageInfo.g",
AbstractHTML := "The package <span class=\"pkgname\">CoReLG</span> contains \
                 functionality for working with real semisimple Lie algebras.",
PackageWWWHome := "http://users.monash.edu/~heikod/corelg/",
Dependencies := rec(
  GAP := ">=4.4",
  NeededOtherPackages:= [ ["sla", ">=0.14"] ],                 
  SuggestedOtherPackages := [ ["GAPDoc", ">= 1.0"] ],
  ExternalConditions := []
),
AvailabilityTest := ReturnTrue,
Autoload := false,

# the banner
BannerString := "CoReLG\n a package for computing with real Lie groups \n by Heiko Dietrich, Paolo Faccin and Willem de Graaf\n",
Keywords := ["real Lie algebras","nilpotent orbits Cartan subalgebras"]
));


