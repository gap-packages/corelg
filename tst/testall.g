LoadPackage("corelg");

# There is a bug in GAP <= 4.10 where factorizing certain polynomials over the
# cyclotomics triggers an incorrect assertion if the assertion level is 2 or
# higher. Unfortunately, the manual tests of corelg trigger the problematic
# case, and START_TEST sets the assertion level to 2 by default. To workaround
# that, we modify START_TEST to only set the assertion level to 1.
#
# Ideally, we'd restrict this hack to GAP < 4.11, by using a check like
#    if not CompareVersionNumbers(GAPInfo.Version, "4.11") then ... fi;
# however, the CI tests we use right now have GAP versions that don't set
# GAPInfo.Version correctly. So for now, we just always enable this workaround.
original_START_TEST := START_TEST;
START_TEST := function( name )
    original_START_TEST( name );
    SetAssertionLevel( 1 );
end;

# run the tests
TestDirectory( DirectoriesPackageLibrary("corelg", "tst"), rec(exitGAP := true,
            testOptions := rec( compareFunction := "uptowhitespace" ) ) );
FORCE_QUIT_GAP(1);
