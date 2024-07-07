#####################################################################################
#
#  PackageInfo.g                    Heiko Dietrich, Paolo Faccin, and Willem de Graaf
#
#
#  The package CoReLG is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by the Free
#  Software Foundation; either version 2 of the License, or (at your option) any
#  later version.

SetPackageInfo( rec(
PackageName := "CoReLG",
Subtitle := "Computing with real Lie algebras",
Version := "1.57",
Date := "07/07/2024", # this is in dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [

  rec(
  LastName := "Dietrich",
  FirstNames := "Heiko",
  IsAuthor := true,
  IsMaintainer := true,
  Email := "heiko.dietrich@monash.edu",
  WWWHome := "http://users.monash.edu.au/~heikod/",
  PostalAddress :=
    """School of Mathematics
    Monash University
    Wellington Road 1
    VIC 3800, Melbourne, Australia""",
  Place := "Melbourne",
  Institution := "School of Mathematical Sciences, Monash University"
  ),

  rec(
  LastName := "Faccin",
  FirstNames := "Paolo",
  IsAuthor := true,
  IsMaintainer := true,
  Email := "paolofaccin86@gmail.com",
  PostalAddress :=
    """Dipartimento di Matematica
    Via Sommarive 14
    I-38050 Povo (Trento), Italy
    """,
  Place := "Trento",
  Institution := "Dipartimento di Matematica, University of Trento"
  ),

  rec(
  LastName := "de Graaf",
  FirstNames := "Willem",
  IsAuthor := true,
  IsMaintainer := true,
  Email := "degraaf@science.unitn.it",
  WWWHome := "https://www.science.unitn.it/~degraaf",
    PostalAddress :=
    """Dipartimento di Matematica
    Via Sommarive 14
    I-38050 Povo (Trento), Italy
    """,
  Place := "Trento",
  Institution := "Dipartimento di Matematica, University of Trento"
  )

],
Status := "accepted",
CommunicatedBy := "Bettina Eick (Braunschweig)",
AcceptDate := "01/2014",

PackageWWWHome  := "https://gap-packages.github.io/corelg/",
README_URL      := Concatenation( ~.PackageWWWHome, "README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "PackageInfo.g" ),
SourceRepository := rec(
    Type := "git",
    URL := "https://github.com/gap-packages/corelg",
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/corelg-", ~.Version ),
ArchiveFormats := ".tar.gz",

PackageDoc := rec( BookName  := "CoReLG" ,
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Computing with real Lie algebras",
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

TestFile := "tst/testall.g",

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
      under grant agreement no 271712, and from the Australian Research
      Council, grantor code DE140100088 and DP190100317.
      """,
    Copyright := "&copyright; 2013-2019 Heiko Dietrich, Paolo Faccin, and Willem de Graaf",
  ),
),

));
