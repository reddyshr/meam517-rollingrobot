#!/bin/perl -w

use strict;

use FindBin;
use lib "$FindBin::Bin/../lib/snadiopt";

use FortranUtils;
use Getopt::Long;
use File::Copy;

use vars qw(%usrfun
            %usrini
            $userfg
            $problem
            $program
            $makefile
            $submakefile
            $refresh_makefile
            $make_defaults
            %sizes
            $VERSION
);
$VERSION = "0.01";

%usrfun = (name => 'usrfun');
%usrini = (name => 'usrini');
$userfg = 'userfg';
$program = 'a.out';
$problem = 'unnamed';
$makefile = 'GNUmakefile';
$refresh_makefile = '';
$make_defaults="$FindBin::Bin/../lib/snadiopt/make_defaults";
%sizes = ();

use vars qw( $opt_help $opt_version $opt_merge);

parse_argv();
read_users_fortran_files(@ARGV);
write_data_contents( \%usrfun );
write_adifor_auxilaries( \%usrfun );
write_admain( \%usrfun );

# End main program

# parse_argv
#
# Get options from the argument list and set the values of variables
# appropriately.

sub parse_argv {
    my $options_ok =
        GetOptions('help'               => \$opt_help,
                   'refresh-makefile'   => \$refresh_makefile,
               'o=s'                => \$program,
               'usrfun=s'           => \$usrfun{name},
               'makefile=s'         => \$makefile,
               'usrini=s'           => \$usrini{name},
               'version'            => \$opt_version,
               'merge'              => \$opt_merge);

    unless ($options_ok) {
        user_help();
        exit( 1 );
    }
    if ($opt_version) {
        print "snadiopt.pl, version $VERSION.\n";
        exit( 0 );
    }
    if ($opt_help) {
        user_help();
        exit( 0 );
    }
    # If the user has given this problem a name, use it.
    unless ( $program eq 'a.out' ) {
        $problem = $program;
    }

    unless( $makefile ) {
        $makefile = "GNUmakefile";
    }
    $submakefile = "${problem}_submake";
}

#   user_help - Print usage information.
sub user_help {
    print &q(<<"    END");
    :Usage: snadiopt.pl [switches] file1.f [file2.f]
    :   -help                 Print this message.
    :   -version              Print the version number of snadiopt.pl.
    :   -o PROGRAM            The optimization problem (and binary executable)
    :                         will be named PROGRAM. (default: a.out)
    :   -makefile MAKEFILE    The output makefile will be named MAKEFILE.
    :                         (default: PROGRAM_submake or unnamed_submake
    :                         if PROGRAM is not specified.)
    :   -refresh-makefile     Create MAKEFILE, even if it already exists.
    :                         Unless given this option, the script will not
    :                         overwrite an existing MAKEFILE.
    :   -usrfun NAME          The FORTRAN subroutine named NAME computes the
    :                         functions needed in this optimization problem.
    :                         (default: usrfun)
    :   -usrini NAME          The FORTRAN subroutine named NAME initializes
    :                         this optimization problem. (default: usrini)
    :   -merge                Merge changes between prob_main.f.orig and
    :                         the current version of prob_main.f into the newly
    :                         generated prob_main.f. Do the same for
    :                         prob_submake.
    END
}

# read_users_fortran_files
#
# Scan the user's supplied FORTRAN files for the sizes of certain arrays
# and the names of the independent and dependent variables.

sub read_users_fortran_files {
    my (@files) = @_;
    my %found;

    my $file;
    foreach $file (@files) {
        open FH, $file
            or die "Can't open $file: $!\n";

        # Have to handle the silly Fortran continuation syntax by
        # doing a one line look-ahead

        my $nextline = <FH>;
        # Skip comment lines
        while ($nextline && $nextline =~ /^[*cC]/) {
            $nextline = <FH>;
        }
        while ($nextline) {
            my $line;
            ($line, $nextline) = freadahead( $nextline, *FH );
            # Allow, but don't preserve, whitespace
            $line =~ s/\s//g;

            if ($line =~ /^subroutine/i) {
                my @namelist = fsubsplit($line);
                my $name = shift @namelist;

                if ($name =~ /^$usrini{name}$/i) {
                    not $found{usrini}++
                    or die "There are two subroutines" .
                        " with name $usrini{name}.\n";

                    $nextline = read_usrini($line, $nextline, *FH, \@namelist);
                    $usrini{file} = $file;
                    $usrini{args} = \@namelist
                } elsif ( $name =~ /^$usrfun{name}$/i ) {
                    not $found{usrfun}++
                    or die "There are two subroutines" .
                        " with name $usrfun{name}.\n";

                    $usrfun{file} = $file;
                    $usrfun{args} = \@namelist;
                }
            }
        }
        close FH;
    } # foreach file

    $found{usrfun}
    or die "The optimization function, '$usrfun{name}', was not defined.\n";
    $found{usrini}
    or die "The data initialization routine, " .
        "'$usrini{name}', was not defined.\n";
}

# write_data_contents
#
# Creates some files for the optimization problem ${problem} based on the
# contents of the data section of this file.
#
# We scan the contents of DATA, and write out each templated file as it is
# found (in other words, the file templates in DATA need not be in any
# particular order).
#
# The single argument to write_data_contents is a reference to the
# structure (hash) containing everything you need to know about the user's
# function definition routine.

sub write_data_contents{
    my ($sub) = @_;

    while ( <DATA> ) {
        if ( /^\s*
            ([A-Za-z_][\w_]*)           # A token
            \s+                         # whitespace
            <<\s*\(\(([^)]+)\)\)\s*$/x  # << ((token))
                ) {
            # This line marks the begining of a section of DATA
            my ($section, $end_token) = ($1, $2);

            SWITCH: for( $section ) {
                /^dad_dispatch_f\b/
                    and
                        write_ad_dispatch($sub, $end_token,
                                          "${problem}_dense_dispatch.f"),
                        last SWITCH;
                /^sad_dispatch_f\b/
                    and
                        write_ad_dispatch($sub, $end_token,
                                          "${problem}_sparse_dispatch.f"),
                        last SWITCH;
                /^snmain_f\b/
                    and
                        write_snopt_main($sub, $end_token ),
                        last SWITCH;
                /^submakefile\b/
                    and write_submakefile($end_token), last SWITCH;
                /^GNUmakefile\b/
                    and ( (! -f $makefile or $refresh_makefile ) and
                        write_GNUmakefile($end_token) ), last SWITCH;
                die 'Unknown case';
            }
        }
    }
}

sub write_GNUmakefile {
    my ($token) = @_;

    open OUTPUT, ">${makefile}"
        or die "Can't open ${makefile} : $!\n";

    LINE: while ( <DATA> ) {
        last if /^\s*$token/;

        if ( /\${make_defaults}/ ) {
            # Insert the platform-specific defaults file into the submakefile.
            open DEFAULTS, $make_defaults
                or die "Can't open $make_defaults : $!\n";

            while(<DEFAULTS>) {
                s|\${LSNADIOPTLIBDIR}|-L$FindBin::Bin/../lib|g;
                print OUTPUT $_;
            }
            close DEFAULTS;
            next LINE;
        }
        # Timestamp the file, if requested.
        if ( /\${now}/ ) {
            my $now = localtime;
            s/\${now}/${now}/g;
        }
    s/\${makefile}/$makefile/g;

        print OUTPUT $_;
    }

    close OUTPUT;
}

#
# write_submakefile
#
# Creates a file named ${submakefile} that contains all the commands
# needed to build a solver for this problem.
#
# The template for this subordinate makefile is contained in the data
# section of this file. When this subroutine is called, the DATA
# filehandle must point to the first line in this template.
#
# If a file named ${makefile} already exists, this routine moves the
# existing file to ${makefile}.bak. The subroutine always creates a file
# named ${makefile}.orig based on the contents of DATA. If $opt_merge is not
# specified, ${makefile} will be the same as ${makefile}.orig. If $opt_merge
# is specified, this routine will attempt to merge with the old ${makefile}.
#
# This subroutine takes a single arguement, a token that marks the end of
# the makefile template. When this token, preceeded only by whitespace,
# is found on a line, processing of the template stops and the subroutine
# returns with the filehandle DATA pointing to the line after the line
# containing the token.

sub write_submakefile {
    my ($token) = @_;

    my ($usrfun, $usrfun_f) = @usrfun{'name', 'file'};
    my $usrini_f;
    if ($usrini{file} ne $usrfun_f) {
        # if the user initialization routine is in a different file
        # from the function definition routine, set $usrini_f to
        # the name of that file
        $usrini_f = $usrini{file}
    } else {
        # otherwise, let $usrini_f be the null string.
        $usrini_f = '';
    }

    if ( -f "${submakefile}" ) {
        # submakefile already exists. Back it up
        rename( "${submakefile}", "${submakefile}.bak" ) or
            die "Can't backup ${submakefile}.";
    }
    open OUTPUT, ">${submakefile}"
        or die "Can't open ${submakefile} : $!\n";

    LINE: while ( <DATA> ) {
        last if /^\s*$token/;
        s/\${usrfun_f}/${usrfun_f}/g;
        s/\${usrini_f}/${usrini_f}/g;
        s/\${usrfun}/${usrfun}/g;
        s/\${program}/${program}/g;
        s/\${problem}/${problem}/g;
        s/\${userfg}/${userfg}/g;

        # Timestamp the file, if requested.
        if ( /\${now}/ ) {
            my $now = localtime;
            s/\${now}/${now}/g;
        }

        print OUTPUT $_;
    }
    close OUTPUT;

    if( $opt_merge && -r "${submakefile}.orig" && -r "${submakefile}.bak" ) {
        # Try to merge with the old makefile, first do the diff
        my $ed_script =
            qx/diff3 -E ${submakefile} ${submakefile}.orig ${submakefile}.bak/;
        # Verify that the diff worked
        $? == 0 or die "Couldn't use diff3 to do the merge";
        # Overwrite the old ${submakefile}.orig
        copy( "${submakefile}", "${submakefile}.orig" )
            or die "Couldn't create ${submakefile}.orig: $!\n";
        # Now use ed to do the merge
        open( EDPIPE, "|ed - ${submakefile}" )
            or die "Couldn't use ed to do the merge: $!\n";
        print EDPIPE $ed_script;
        print EDPIPE "wq\n";
        close( EDPIPE )
            or die "Couldn't complete the merge of ${submakefile}.\n"
    } else {
        # The user has not requested a merge. Just copy the newly-created
        # submakefile into ${submakefile}.orig
        copy( "${submakefile}", "${submakefile}.orig" )
            or die "Couldn't create ${submakefile}.orig: $!\n";
    }
}


sub write_ad_dispatch {
    my ($sub, $token, $filename) = @_;
    my $usrfun = $sub->{name};

    open OUTFILE, ">$filename"
        or die "Can't open $filename: $!\n";

    while( <DATA> ) {
        last if /^\s*$token/;
        if ( s/\${userfg}/$userfg/g || s/\${usrfun}/$usrfun/g ) {
            $_ = fortwrapline($_);
        }
        print OUTFILE;
    }
    close OUTFILE;
}

# write_snopt_main
#
sub write_snopt_main{
    my ($sub, $token) = @_;

    my ($nName, $nFname, $n, $mCon,
        $lencu, $leniu, $lenru,
        $lencw, $leniw, $lenrw)
        = @sizes{qw(nName nFname n mCon
                    lencu leniu lenru
                    lencw leniw lenrw)};

    my ($call_usrini);

    my ($mtrue, $ntrue) = ('neF', 'n');

    $call_usrini = &q(<<"    END");
    :call $usrini{name}( ObjAdd, ObjRow, Prob,
    :     &     x, bl(1), bu(1), xstate, Names,
    :     &     Fmul, bl(n+1), bu(n+1), Fstate, Fnames,
    :     &     iSpecs, iPrint, iSumm, iErr,
    :     &     cu, iu, ru, cw, iw, rw )
    END


    my $parameters = &q(<<"    END");
    :      parameter( lencw = $lencw, leniw = $leniw, lenrw = $lenrw )
    :      parameter( lencu = $lencu, leniu = $leniu, lenru = $lenru )
    :      parameter( n = $n, neF = $mCon )
    :      parameter( nName = $nName, nFname = $nFname )
    END

    if ( -f "${problem}_main.f" ) {
        rename( "${problem}_main.f", "${problem}_main.f.bak" ) or
            die "Can't backup ${problem}_main.f\n";
    }

    open OUTFILE, ">${problem}_main.f"
        or die "Can't open ${problem}_main.f: $!\n";

    while( <DATA> ) {
        last if /^\s*$token/;
        if ( s/\${call_usrini}/${call_usrini}/g ||
            s/\${parameters}/${parameters}/g ) {
                $_ = fortwrap($_);
        } elsif( /\${decl_mtrue}/ ) {
            if ( $mtrue ne 'neF' ) {
                $_ =
                fortwrap("      integer             ${mtrue}, ${ntrue}\n");
            } else {
                $_ = '';
            }
        } elsif( /\${problem}/ ) {
                my $name8 = sprintf "%-8.8s", $program;
                s/\${problem}/$name8/g;
        } else {
            my $needs_wrap = 0;
            $needs_wrap = 1 if s/\${userfg}/${userfg}/g;
            $needs_wrap = 1 if s/\${mtrue}/${mtrue}/g;
            $needs_wrap = 1 if s/\${ntrue}/${ntrue}/g;

            $_ = fortwrap($_) if $needs_wrap;
        }

        print OUTFILE;
    }

    close OUTFILE;

    if( $opt_merge && -r "${problem}_main.f.orig" &&
       -r "${problem}_main.f.bak" ) {
        # We're going to try a merge, first do the diff
        my $ed_script =
   qx/diff3 -E ${problem}_main.f ${problem}_main.f.orig ${problem}_main.f.bak/;
        # Check the status from the diff
        $? == 0 or die "Couldn't use diff3 to do the merge: $!\n";
        # nuke the old version of ${problem}_main.f.orig
        copy( "${problem}_main.f", "${problem}_main.f.orig" )
            or die "Can't create ${problem}_main.f.orig: $!\n";
        # merge with ${problem}_main.f
        open(EDPIPE, '|ed - ${problem}_main.f' )
            or die "Couldn't use ed to do the merge: $!\n";
        print EDPIPE $ed_script;
        print EDPIPE "wq\n";
        close( EDPIPE ) or "The merge of ${problem}_main.f failed: $!\n";

    } else {
        # We aren't doing a merge
        copy( "${problem}_main.f", "${problem}_main.f.orig" )
            or die "Can't create ${problem}_main.f.orig: $!\n";
    }
}

sub write_adifor_auxilaries {
    my ($sub) = @_;

    my $usrfun = $sub->{name};
    my $x = $sub->{args}[4];
    my $f = $sub->{args}[5];

    open ADFFILE, ">${problem}.adf"
        or die "Can't open ${problem}.adf: $!\n";

    print ADFFILE &q(<<"    END");
    :AD_PROG = ${problem}.cmp
    :AD_TOP  = $usrfun
    :AD_IVARS = $x
    :AD_DVARS = $f
    :AD_PMAX = $sizes{n}
    :AD_OUTPUT_DIR = .
    END

    close ADFFILE;
}
sub write_admain {
    my ($sub) = @_;

    my $file  = $sub->{file};
    my $usrfun = $sub->{name};

    my ($nState, $mode, $mCon, $n, $x, $f,
         $cu, $lencu, $iu, $leniu, $ru, $lenru,
         $cw, $lencw, $iw, $leniw, $rw, $lenrw) =
             @{$sub->{args}};

    my $outname = "${problem}_admain.$fortran_extension";

    open OUT, ">$outname" or die "Can't open file $outname: $!.\n";

    print OUT fortwrap(&q(<<"    END"));
    :      program            ${usrfun}_main
    :      implicit           none
    :      integer            $nState, $mode, $mCon, $n
    :      parameter ( $nState = 1, $mode = 1, $mCon = 1, $n = 1 )
    :      double precision   $x($n), $f($mCon)
    :
    :      integer            $lencu, $leniu, $lenru, $lencw, $leniw, $lenrw
    :      parameter ( $lencu = 1, $leniu = 1, $lenru = 1,
    :     &            $lencw = 1, $leniw = 1, $lenrw = 1 )
    :
    :      character*8        $cu($lencu), $cw($lencw)
    :      integer            $iu($leniu), $iw($leniw)
    :      double precision   $ru($lenru), $rw($lenrw)
    :
    :      call  $usrfun(  $nState, $mode,
    :     &              $mCon, $n, $x, $f,
    :     &              $cu, $lencu, $iu, $leniu, $ru, $lenru,
    :     &              $cw, $lencw, $iw, $leniw, $rw, $lenrw )
    :
    :
    :      end
    :
    :      subroutine snset ( buffer, iPrint, iSumm, INFO,
    :     &     cw, lencw, iw, leniw, rw, lenrw )
    :      implicit
    :     &     none
    :      character*(*)
    :     &     buffer
    :      integer
    :     &     iPrint, iSumm, INFO, lencw, leniw, lenrw, iw(leniw)
    :      double precision
    :     &     rw(lenrw)
    :      character*8
    :     &     cw(lencw)
    :      end
    :
    :      subroutine snseti( buffer, ivalue, iPrint, iSumm, INFO,
    :     &     cw, lencw, iw, leniw, rw, lenrw )
    :      implicit
    :     &     none
    :      character*(*)
    :     &     buffer
    :      integer
    :     &     ivalue, iPrint, iSumm, INFO, lencw, leniw, lenrw,
    :     &     iw(leniw)
    :      double precision
    :     &     rw(lenrw)
    :      character*8
    :     &     cw(lencw)
    :      end
    :
    :      subroutine snsetr( buffer, rvalue, iPrint, iSumm, INFO,
    :     &     cw, lencw, iw, leniw, rw, lenrw )
    :      implicit
    :     &     none
    :      character*(*)
    :     &     buffer
    :      integer
    :     &     iPrint, iSumm, INFO, lencw, leniw, lenrw, iw(leniw)
    :      double precision
    :     &     rvalue, rw(lenrw)
    :      character*8
    :     &     cw(lencw)
    :      end
    END
    close OUT;
}

sub q{
    my $string  = $_[0];
    $string =~ s/^\s*://gm;
    #$string =~ s{ (.*)\*/\s*$ }{ sprintf "%-73s*/\n", $1 }gmex;

    return $string;
}
sub read_usrini {
    my ($line, $nextline, $FH, $namelist) = @_;

    scalar(@{$namelist}) == 23
        or die
            "The user inialization routine has the wrong number " .
            "of parameters.\n" .
                "Expected 23 got " . scalar(@{$namelist}) . ".\n";
    my %dim = (
               $namelist->[3]   => 'n',
               $namelist->[7]   => 'nName',
           $namelist->[8]   => 'mCon',
               $namelist->[12]  => 'nFname',
           $namelist->[17]  => 'lencu',
           $namelist->[18]  => 'leniu',
           $namelist->[19]  => 'lenru',
           $namelist->[20]  => 'lencw',
           $namelist->[21]  => 'leniw',
           $namelist->[22]  => 'lenrw');

    @sizes{values(%dim)} = ();

    my %params = ();

    LINE: while( $nextline ) {
        ($line, $nextline) = freadahead( $nextline, $FH );
        # Allow, but don't preserve, whitespace
        $line =~ s/\s//g;
        last LINE if $line =~ /^end$/i;

        if ($line =~ /^parameter/i) {
            my @newParam = fparamsplit( $line );
            my $i;
            for ( $i = 0; $i < scalar(@newParam); $i += 2 ) {
                $params{$newParam[$i]} = $newParam[$i+1];
            }
        } elsif ( $line =~
                 /^integer |
                  ^doubleprecision |
                  ^character |
                  ^dimension/ix ) {
            my @decl = fdeclsplit($line);
            shift @decl; #ignore the data type (perhaps) check this some day

            my $decl;
            DECL: foreach $decl (@decl) {
                $decl =~ /^($fortran_token)(?:\((.*)\))?/
                    or die 'Syntax error in declaration';

                if ( $dim{$1} ) {
                    my ($name, $arg) = ($1, $2);

                    unless($arg) {
                        # There should be an argument giving the array
                        # dimension. There isn't, but this may not be an
                        # error. They may use a dimension statement later
                        next DECL;
                    }
                    while( $arg ) {
                        if ($arg =~ /^\+?\d+$/) {
                            # $arg is an integer
                            $sizes{$dim{$name}} = $arg;
                            $arg = '';
                        } elsif ($arg =~ /^$fortran_token$/) {
                            # $arg may be a parameter
                            exists $params{$arg}
                            or die "$arg is not a known parameter\n";

                            $arg = $params{$arg};
                        } else {
                            die "$arg is not a number or a parameter name\n";
                        }
                    }
                }
            }
        }
    }

    foreach (keys %sizes) {
        unless( $sizes{$_} ) {
            die "$_ was not defined\n";
        }
    }
    return $nextline;
}

__DATA__

dad_dispatch_f << ((END_DAD_DISPATCH))
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     dad_dispatch.f (automatically generated by snadiopt.pl)
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nwadmx
     &   ( imtx, limtx, isize, morei )

      implicit
     &     none
      integer
     &     limtx, imtx(limtx), isize, morei

      call nwddmx
     &   ( imtx, limtx, isize, morei )

      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userfg
     &   ( Status,
     &     nv, v,
     &     needF, neF, F,
     &     needG, neG, G,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit           none
      integer            Status, needF, needG
      integer            nv, neF
      double precision   v(nv), F(neF)

      integer            neG
      double precision   G(neG)

      integer            lencu, leniu, lenru
      integer            lencw, leniw, lenrw

      character*8        cu(lencu), cw(lencw)
      integer            iu(leniu), iw(leniw)
      double precision   ru(lenru), rw(lenrw)

      integer            limtx, lrmtx

      external g_${usrfun}
      integer             iiadmx, liadmx, iradmx, lradmx
      integer             iiwork, liwork, irwork, lrwork

      integer             nsnr
      parameter         ( nsnr = 500 )

      call getnlp( limtx, lrmtx,
     &     iiadmx, liadmx, iradmx, lradmx,
     &     iiwork, liwork, irwork, lrwork,
     &     iw(nsnr + 1), leniw - nsnr )

      call ddfnfg( Status,
     &     nv, v,
     &     needF, neF, F,
     &     needG, neG, G,
     &     g_${usrfun},
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     iw(nsnr + iiadmx), liadmx,
     &     rw(nsnr + iradmx), lradmx )

      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine iniad
     &   ( v0, nv, blv, buv,
     &     bl, neF, bu,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     imtx, limtx, rmtx, lrmtx,
     &     lusdi, lusdr, morei, morer )

      implicit           none

      integer            nv
      double precision   v0(nv), blv(nv), buv(nv)
      integer            neF
      double precision   bl(neF), bu(neF)
      integer            ladtlc
      external           g_${usrfun}, ladtlc, adtlc
      integer            lencu, lencw
      integer            leniu, leniw, limtx
      integer            lenru, lenrw, lrmtx
      character*8        cu(lencu), cw(lencw)
      integer            iu(leniu), iw(leniw), imtx(limtx)
      double precision   ru(lenru), rw(lenrw), rmtx(lrmtx)
      integer            lusdi, lusdr, morei, morer

      call inidad
     &   ( v0, nv, blv, buv,
     &     bl, neF, bu,
     &     g_${usrfun}, adtlc, ladtlc,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     imtx, limtx, rmtx, lrmtx,
     &     lusdi, lusdr, morei, morer )

      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      function ladtlc()
      integer  ladtlc

      ladtlc = 0
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine adtlc ( iAfun, neA, jAvar, A )

      implicit           none
      integer            neA, iAfun(neA), jAvar(neA)
      double precision   A(neA)
*     Relax
      end

END_DAD_DISPATCH

sad_dispatch_f << ((END_SAD_DISPATCH))
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     sad_dispatch.f (automatically generated by snadiopt.pl)
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nwadmx( imtx, limtx, isize, morei )
      implicit       none
      integer        limtx, imtx(limtx), isize, morei

      call nwsdmx( imtx, limtx, isize, morei )

      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine userfg
     &   ( Status,
     &     nv, v,
     &     needF, neF, F,
     &     needG, neG, G,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw)

      implicit           none

      integer            Status, needF, needG
      integer            nv, neF
      double precision   v(nv), F(neF)

      integer            neG
      double precision   G(neG)

      integer            lencu, leniu, lenru
      integer            lencw, leniw, lenrw

      character*8        cu(lencu), cw(lencw)
      integer            iu(leniu), iw(leniw)
      double precision   ru(lenru), rw(lenrw)

      integer            limtx, lrmtx

      external g_${usrfun}
      integer            iiadmx, liadmx, iradmx, lradmx
      integer            iiwork, liwork, irwork, lrwork

      integer            nsnr
      parameter        ( nsnr = 500 )


      call getnlp
     &   ( limtx, lrmtx,
     &     iiadmx, liadmx, iradmx, lradmx,
     &     iiwork, liwork, irwork, lrwork,
     &     iw(nsnr + 1), leniw - nsnr )

      call sdfnfg
     &   ( Status,
     &     nv, v,
     &     needF, neF, F,
     &     needG, neG, G,
     &     iw(nsnr + iiwork), rw(nsnr + irwork),
     &     g_${usrfun},
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     iw(nsnr + iiadmx), liadmx,
     &     rw(nsnr + iradmx), lradmx )

      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine iniad
     &   ( v0, nv, blv, buv,
     &     bl, neF, bu,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     imtx, limtx, rmtx, lrmtx,
     &     lusdi, lusdr, morei, morer )

      implicit           none

      integer            nv
      double precision   v0(nv), blv(nv), buv(nv)
      integer            neF
      double precision   bl(neF), bu(neF)
      integer            ladtlc
      external           g_${usrfun}, ladtlc, adtlc
      integer            lencu, lencw
      integer            leniu, leniw, limtx
      integer            lenru, lenrw, lrmtx
      character*8        cu(lencu), cw(lencw)
      integer            iu(leniu), iw(leniw), imtx(limtx)
      double precision   ru(lenru), rw(lenrw), rmtx(lrmtx)
      integer            lusdi, lusdr, morei, morer

      call inisad
     &   ( v0, nv, blv, buv,
     &     bl, neF, bu,
     &     g_${usrfun}, adtlc, ladtlc,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     imtx, limtx, rmtx, lrmtx,
     &     lusdi, lusdr, morei, morer )

      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      function ladtlc( )
      integer  ladtlc

      ladtlc = 0
      end

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine adtlc ( iAfun, neA, jAvar, A )

      implicit           none
      integer            neA, iAfun(neA), jAvar(neA)
      double precision   A(neA)
*     Relax
      end

END_SAD_DISPATCH

snmain_f << ((END_SNMAIN))
      program            snmain

      implicit           none

      external           ${userfg}

      integer            n, neF, nName, nFname
      integer            lencw, leniw, lenrw, lencu, leniu, lenru

*
${parameters}
*

      character*8        Names(nName), Fnames( nFname )

      integer            nsnr
      parameter        ( nsnr = 500 )

      character*8        cu(lencu),          cw(lencw)
      integer            iu(leniu),          iw(leniw)
      double precision   ru(lenru),          rw(lenrw)

      integer            nbnd
      parameter        ( nbnd = n + neF )

      double precision   bl(nbnd), bu(nbnd)
      integer            xstate(n)
      double precision   x(n), xmul(n)
      integer            Fstate(neF)
      double precision   F(neF), Fmul(neF)

      double precision   nInfty,           pInfty
      parameter         (nInfty = -1.0d20, pInfty = 1.0d20)
      double precision   zero
      parameter         (zero = 0.0d0)

      character*8        Prob

      integer            INFO, mincw, miniw, minrw
      integer            nS, nInf
      double precision   sInf, ObjAdd

      integer            iSpecs, iPrint, iSumm, iErr

      integer            iiAfun, ijAvar, lenA, neA, iA
      integer            iiGfun, ijGvar, lenG, neG, mm

      integer            lPrint, lSumm
      integer            lusdi, lusdr, moreiw, morerw
      integer            ObjRow, iCold
      ${decl_mtrue}

      integer            iiadmx, liadmx, iradmx, lradmx
      integer            linlp,  lrnlp
      integer            iiwork, liwork, irwork, lrwork
      integer            lenia,  lenra
*
*     Data statements
*
      data bl          / nbnd * nInfty /
      data bu          / nbnd * pInfty /
      data xstate      / n * 0 /
      data x           / n * zero /
      data xmul        / n * zero /
      data Fstate      / neF * 0 /
      data F           / neF * zero/
      data Fmul        / neF * zero/

      data ObjAdd      / zero /
      data ObjRow      / 1 /
      data neA, neG    / 0, 0 /
      data Prob        /'${problem}'/

      data iSpecs, iPrint, iSumm, iErr /0, 0, 6, 0/
      data         lPrint, lSumm       /   0, 0   /
*
*     End Data statements
*

*     ------------------------------------------------------------------
*     First,  snInit MUST be called to initialize optional parameters
*     to their default values.
*     ------------------------------------------------------------------
      call snInit
     &  ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

      ${call_usrini}

*     ------------------------------------------------------------------
*     Reset iPrint and iSumm from usrini.
*     ------------------------------------------------------------------
      call snseti
     &   ( 'iw 12 = ', iPrint, lPrint, lSumm, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call snseti
     &   ( 'iw 13 = ', iSumm , lPrint, lSumm, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

      lenia          = leniw - nsnr
      lenra          = lenrw - nsnr
      lusdi          = 0
      lusdr          = 0

      call newnlp( iw(nsnr + 1), leniw - nsnr, lusdi, moreiw )
      if ( moreiw .ne. 0 ) goto 600

      call ininlp( x, ${ntrue}, bl(1), bu(1),
     &     bl(n+1), ${mtrue}, bu(n+1),
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw,
     &     iw(nsnr + 1), lenia, rw(nsnr + 1), lenra,
     &     lusdi, lusdr, moreiw, morerw )

      if (moreiw .ne. 0  .or.  morerw .ne. 0) goto 600
      lenia = lusdi
      lenra = lusdr

*     Get the copy of the Jacobian (it is a copy because snopt
*     scrambles it.

      call getnlp( linlp, lrnlp,
     &     iiadmx, liadmx, iradmx, lradmx,
     &     iiwork, liwork, irwork, lrwork,
     &     iw(nsnr + 1), lenia )

      call getmtx( iiAfun, neA, ijAvar, iA,
     &     iiGfun, neG, ijGvar, mm,
     &     iw(nsnr + iiadmx), liadmx, rw(nsnr + iradmx), lradmx )

*     Adjust iiAfun, etc. to be indices into iw and rw
      iiAfun = iiAfun + nsnr + iiadmx - 1
      ijAvar = ijAvar + nsnr + iiadmx - 1
      iiGfun = iiGfun + nsnr + iiadmx - 1
      ijGvar = ijGvar + nsnr + iiadmx - 1

      iA     = iA     + nsnr + iradmx - 1

      call snseti
     &   ( 'User real    workspace', lusdr + nsnr,
     &     lPrint, lSumm, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call snseti
     &   ( 'User integer workspace', lusdi + nsnr,
     &     lPrint, lSumm, INFO,
     &     cw, lencw, iw, leniw, rw, lenrw )

*     ------------------------------------------------------------------
*     Read a Specs file (Optional).
*     ------------------------------------------------------------------
      if (iSpecs .gt. 0) then
         call snSpec
     &     ( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

         if (INFO .ne. 101) then
             write(iErr, *) 'Trouble reading the Specs file'
             stop 1
         end if
      end if

      iCold = 0
      lenA  = max( neA, 1 )
      lenG  = max( neG, 1 )

      call snOptA
     &   ( iCold, ${mtrue}, ${ntrue}, nName, nFname,
     &     ObjAdd, ObjRow, Prob, ${userfg},
     &     iw(iiAfun), iw(ijAvar), lenA, neA, rw(iA),
     &     iw(iiGfun), iw(ijGvar), lenG, neG,
     &     bl(1), bu(1), Names, bl(n+1), bu(n+1), FNames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (INFO .ge. 10) goto 900

      stop

*     Begin error handling

 600  INFO  = 1
      miniw = leniw + moreiw
      minrw = lenrw + morerw

 900  if ( INFO .eq. 1 ) then

*        Memory allocation errors

         write (iErr, *) 'There is not enough storage ',
     &        'to start solving the problem'

         write (iErr, *) 'Total integer workspace should be ',
     &        'significantly more than ', miniw

         write (iErr, *) 'Total real    workspace should be ',
     &        'significantly more than ', minrw

         stop 1
      else
         write (iErr, *) Prob, ' terminated with INFO = ',
     &        INFO
         stop 2
      end if

      end

END_SNMAIN

GNUmakefile << ((END_GNU_MAKEFILE))
# ${makefile}
#
# Created:         ${now} by snadiopt.pl
# Current Version: ${now}


${make_defaults}

# Set environment variables used by the ADIFOR pre-processor
export AD_HOME AD_OS AD_LIB PATH

# ADIFORLIBS - Libraries that must be linked to ADIFOR-generated code
ADIFORLIBS = \
    $(AD_LIB)/lib/ReqADIntrinsics-$(AD_OS).o \
    $(AD_LIB)/lib/libADIntrinsics-$(AD_OS).a

SPARSE_LIB     = $(AD_LIB)/lib/libSparsLinC-$(AD_OS).a

AD_FLAVOR := sparse

# FTNCHEK - The name of a program to use to verify the consistency of
#    the fortran files. Unsurprisingly, this program is usually named
#    ftnchek. The ftnchek program is not distributed with snopt,
#    but is readily available (and highly recommended). Try looking at
#    ftp://netlib.org/fortran/ftnchek.tar.gz.
#
#    'make check' will perform the consistency check.
FTNCHEK = ftnchek

# FTNCHEKFLAGS - Flags used by the ftnchek program. The default flags
#    are valid for ftncheck v2.09, and turn off unavoidable warnings
#    in snopt-adifor generated programs.
FTNCHEKFLAGS = -noextern -pretty=all -pretty=no-embedded-space \
    -quiet -array=0 -usage=331

# Make all the first target
all:

include $(wildcard *_submake)

# make check - checks the consistency of the fortran files used to
#     build each program. Requires that you have ftnchek installed.
check:

# make clean - removes object files and some common "garbage" files
#   as "core" files.
clean:

# veryclean - removes more files generated by the build process,
#             including the executable. Removes the output file.
#
veryclean: clean adifor-veryclean

# make distclean - cleans up for distribution.
#     The adifor-generated fortran files are considered part of the
#     distribution.
distclean: clean adifor-clean snadiopt-clean

# make maintainer-clean - deletes everything that can be rebuilt
maintainer-clean: veryclean adifor-veryclean snadiopt-veryclean
	@echo 'This command is intended for maintainers to use; it'
	@echo 'deletes files that may need special tools to rebuild.'
	-rm GNUmakefile

# make adifor-clean - Removes ADIFOR auxiliarly files
adifor-clean:
	-rm -rf AD_cache

# make adifor-veryclean - Removes ADIFOR auxiliarly files and
#     the generated fortran files  (and *.cmp).
adifor-veryclean: adifor-clean

# snadiopt-clean - Removes auxiliary files generated by snadiopt
#   for use with all the programs
snadiopt-clean:

# snadiopt-veryclean - Removes all files generated by snadiopt.pl
#   for use with all the programs
snadiopt-veryclean: snadiopt-clean

# Cancel the rules for building an executable from a .o or .f file.
# This helps protect against typos.

% : %.o
% : %.f

.PHONY : check
.PHONY : clean veryclean
.PHONY : distclean maintainer-clean
.PHONY : adifor-clean adifor-veryclean
.PHONY : snadiopt-clean snadiopt-veryclean
END_GNU_MAKEFILE

submakefile << ((END_SUBMAKEFILE))
all: ${program}
${problem}-all: ${program}

# ${problem}_USER_LIBS - Additional libraries that a user might want
#                        to supply.
${problem}_USER_LIBS =

# Sparse Adifor
${problem}_SPARSE_LIBS  =  \
    $(SPARSESNADIOPTLIBS) $(ADIFORLIBS) $(SPARSE_LIB) $(BLAS)
${problem}_SPARSE_SOURCE  = ${problem}_main.f ${problem}_sparse_dispatch.f  \
    $(${problem}_AD_SOURCE) $(${problem}_AD_G_SOURCE)
${problem}_SPARSE_OBJECTS = $(${problem}_SPARSE_SOURCE:.f=.o)

# Dense Adifor
${problem}_DENSE_LIBS  =  \
    $(DENSESNADIOPTLIBS) $(ADIFORLIBS) $(BLAS)
${problem}_DENSE_SOURCE  = ${problem}_main.f ${problem}_dense_dispatch.f \
    $(${problem}_AD_SOURCE) $(${problem}_AD_G_SOURCE)
${problem}_DENSE_OBJECTS = $(${problem}_DENSE_SOURCE:.f=.o)

# ${problem}_LIBS    - All libraries which must be linked to the final program.
# ${problem}_SOURCE  - All source files for ${problem}
# ${problem}_OBJECTS - All object files for ${problem}
ifeq ($(AD_FLAVOR),sparse)

${problem}_LIBS    = $(${problem}_USER_LIBS) $(${problem}_SPARSE_LIBS)
${problem}_SOURCE  = $(${problem}_SPARSE_SOURCE)
${problem}_OBJECTS = $(${problem}_SPARSE_OBJECTS)

else

${problem}_LIBS    = $(${problem}_USER_LIBS) $(${problem}_DENSE_LIBS)
${problem}_SOURCE  = $(${problem}_DENSE_SOURCE)
${problem}_OBJECTS = $(${problem}_DENSE_OBJECTS)

endif

# OBJECTS - object files for ${problem}
${problem}_OBJECTS = $(${problem}_SOURCE:.f=.o)

# AD_SOURCE - Fortran files that will be differentiated by ADIFOR
${problem}_AD_SOURCE = ${usrfun_f}

# AD_G_SOURCE - Fortran files produced by the ADIFOR pre-processor
${problem}_AD_G_SOURCE = $(${problem}_AD_SOURCE:%.f=g_%.f)


# AD_OTHER_FILES - Fortran files that need to be passed to ADIFOR to make
#     a complete program. These files differ from the files in
#     AD_SOURCE in that for each filename.f in AD_SOURCE, ADIFOR will
#     generate a correspoing file named g_filename.f which contains
#     source for ${program}. AD_OTHER_FILES should only contain source
#     files for which ADIFOR will not generate output files with names
#     like g_filename.f.
${problem}_AD_OTHER_FILES = ${problem}_admain.f

# AD_AUX - Contains the names of the copious auxilary files which
# adifor generates.
${problem}_AD_AUX     = \
    $(${problem}_AD_SOURCE:%.f=.%.f) \
    $(${problem}_AD_OTHER_FILES:%.f=.%.f) \
    $(${problem}_AD_G_SOURCE:%.f=%.A) \
    $(${problem}_AD_G_SOURCE:%.f=%.aux)

# Build the program
${program} : $(${problem}_OBJECTS)
	$(LINK.f) -o $@ $(${problem}_OBJECTS) $(${problem}_LIBS) \
	    $(LDLIBES) $(LOADLIBES)

$(${problem}_AD_G_SOURCE) : $(${problem}_AD_OTHER_FILES) \
	$(${problem}_AD_SOURCE) ${problem}.adf ${problem}.cmp
	$(ADIFOR) AD_SCRIPT=${problem}.adf AD_FLAVOR=$(AD_FLAVOR)
	touch $(${problem}_AD_G_SOURCE)
	# Make all files up-to-date, even if they weren't
	# changed by ADIFOR

# make ${problem}.cmp -
#     Generates ${problem}.cmp from the variables internal to this makefile.
${problem}.cmp : $(MAKEFILE)
	ls -1 $(${problem}_AD_SOURCE) $(${problem}_AD_OTHER_FILES) \
	    > ${problem}.cmp

# make ${problem}-check - checks the consistency of the fortran files used to
#     build ${problem}. Requires that you have ftnchek installed.
${problem}-check:
	$(FTNCHEK) $(FTNCHEKFLAGS) $(${problem}_SOURCE)

# make ${problem}-clean - removes object files
#     and some common "garbage" files such as "core" files.
${problem}-clean: ${problem}-adifor-clean
	-rm -f $(${problem}_DENSE_OBJECTS) $(${problem}_SPARSE_OBJECTS)
	-rm -f *~ core

# make ${problem}-adifor-clean - removes ADIFOR auxiliarly files
${problem}-adifor-clean:
	-rm -f $(${problem}_AD_AUX)

# make ${problem}-adifor-veryclean - removes ADIFOR auxiliarly files and
#     the generated fortran files  (and ${problem}.cmp).
${problem}-adifor-veryclean: ${problem}-adifor-clean
	-rm -f $(${problem}_AD_G_SOURCE) ${problem}.cmp

# ${problem}-snadiopt-clean - removes auxiliary files generated by snadiopt
#   for use with ${problem}
${problem}-snadiopt-clean:
	-rm -f ${problem}_submake.orig ${problem}_submake.bak
	-rm -f ${problem}_main.f.orig        ${problem}_main.f.bak

# ${problem}-snadiopt-veryclean - removes all files generated by snadiopt.pl
#   for use with ${problem}
${problem}-snadiopt-veryclean: ${problem}-snadiopt-clean
	-rm -f ${problem}_admain.f ${problem}.adf ${problem}.cmp \
	    ${problem}_dense_dispatch.f ${problem}_sparse_dispatch.f \
	    ${problem}_main.f ${problem}_submake

# make ${problem}-veryclean - removes more files generated by the build
#   process, including the executable. Removes the output file.
${problem}-veryclean: ${problem}-clean ${problem}-adifor-veryclean
	-rm -f ${program}
	-rm -f ${program}.out

# make ${problem}-distclean - cleans up for distribution.
#     The adifor-generated fortran files are considered part of the
#     distribution.
${problem}-distclean: ${problem}-clean ${problem}-adifor-clean \
    ${problem}-snadiopt-clean
	-rm -f *~ core
	-rm -f ${program}
	-rm -f ${program}.out

# make ${problem}-maintainer-clean - deletes everything that can be rebuilt
${problem}-maintainer-clean: ${problem}-veryclean ${problem}-adifor-veryclean \
        ${problem}-snadiopt-veryclean

check: ${problem}-check
clean:  ${problem}-clean
veryclean: ${problem}-veryclean
distclean: ${problem}-distclean
maintainer-clean: ${problem}-maintainer-clean

adifor-clean: ${problem}-adifor-clean
adifor-veryclean: ${problem}-adifor-veryclean

snadiopt-clean: ${problem}-snadiopt-clean
snadiopt-veryclean: ${problem}-snadiopt-veryclean

.PHONY : ${problem}-check
.PHONY : ${problem}-clean ${problem}-veryclean
.PHONY : ${problem}-distclean ${problem}-maintainer-clean
.PHONY : ${problem}-adifor-clean ${problem}-adifor-veryclean
.PHONY : ${problem}-snadiopt-clean ${problem}-snadiopt-veryclean

END_SUBMAKEFILE

=head1 NAME

snadiopt.pl - An interface between the numerical
optimization package B<snopt> and the automatic differentiation
package B<ADIFOR>.


=head1 SYNOPSIS

B<snadiopt.pl>
S<[ B<-help> ]> S<[ B<-version> ]> S<[ B<-o> I<PROGRAM>] >
S<[ B<-makefile> I<MAKEFILE> ]> S<[ B<-refresh-makefile> ] >
S<[ B<-usrfun> I<NAME> ]> S<[ B<-usrini> I<NAME> ]> S<[ B<-merge>]>
file1.f S<[ file2.f ] >

=head1 DESCRIPTION

=head2 Switches

The command line options are:

=over 4

=item B<-help>

Print a summary of command line options.

=item B<-version>

Print the version number of snadiopt.pl.

=item B<-o> I<PROGRAM>

The optimization problem (and binary executable) will be named PROGRAM.
(default: a.out)

=item B<-makefile> I<MAKEFILE>

The output makefile will be named MAKEFILE. (default: GNUmakefile)

=item B<-refresh-makefile>

Create MAKEFILE, even if it already exists.  Unless given this option,
the script will not overwrite an existing MAKEFILE.

=item B<-usrfun> I<NAME>

The FORTRAN subroutine named NAME computes the functions needed in this
optimization problem.  (default: usrfun)

=item B<-usrini> I<NAME>

The FORTRAN subroutine named NAME initializes this optimization problem.
(default: usrini)

=item B<-merge>

Merge changes between snmain.f.orig and the current version of
snmain.f into the newly generated snmain.f. Do the same for the
makefile.

=back

=cut
