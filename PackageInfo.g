#####################################################################################
#
# PackageInfo.g Heiko Dietrich, Paolo Faccin, and Willem de Graaf
#
#
# The package CoReLG is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software Foundation;
# either version 2 of the License, or (at your option) any later version.

# Details may have to be corrected...


SetPackageInfo( rec(
PackageName := "CoReLG",
Subtitle := "Computation with real Lie groups",
Version := "1.51",
Date := "15/07/2019", # this is in dd/mm/yyyy format
License := "GPL-2.0-or-later",
                
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

PackageWWWHome := "https://gap-packages.github.io/corelg/",
README_URL := Concatenation( ~.PackageWWWHome, "README.md" ),
PackageInfoURL := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
SourceRepository := rec(
Type := "git",
URL := "https://github.com/gap-packages/corelg",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
ArchiveURL := Concatenation( ~.SourceRepository.URL,
"/releases/download/v", ~.Version,
"/corelg-", ~.Version ),
ArchiveFormats := ".tar.gz",

PackageDoc := rec( BookName := "CoReLG" ,
ArchiveURLSubset := ["doc"],
HTMLStart := "doc/chap0.html",
PDFFile := "doc/manual.pdf",
SixFile := "doc/manual.six",
LongTitle := "Computing with real Lie groups",
Autoload := false
),

AbstractHTML := "The package <span class=\"pkgname\">CoReLG</span> contains \
functionality for working with real semisimple Lie algebras.",

Dependencies := rec(
GAP := ">=4.8",
NeededOtherPackages:= [ ["sla", ">=1.5"] ],
SuggestedOtherPackages := [ ],
ExternalConditions := []
),
AvailabilityTest := ReturnTrue,
Autoload := false,

# the banner
BannerString := "CoReLG\n a package for computing with real Lie groups \n by Heiko Dietrich, Paolo Faccin and Willem de Graaf\n",
Keywords := ["real Lie algebras","nilpotent orbits Cartan subalgebras"],

AutoDoc := rec(
TitlePage := rec(
Version := Concatenation( "Version ", ~.Version ),
Abstract := """
This package provides functions for computing with various
aspects of the theory of real simple Lie algebras.
""",
Acknowledgements := """
The research leading to this package has received funding from
the European Union's Seventh Framework Program FP7/2007-2013
under grant agreement no 271712.
""",
Copyright := "&copyright; 2014 Heiko Dietrich, Paolo Faccin, and Willem de Graaf",
),
),

));

