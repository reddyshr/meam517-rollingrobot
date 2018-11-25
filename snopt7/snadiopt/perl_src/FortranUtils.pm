package FortranUtils;

use strict;
use vars qw($VERSION @ISA @EXPORT);

use vars qw($fortran_extension $fortran_token $fortran_number);

require Exporter;

@ISA = qw(Exporter);
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.
@EXPORT = qw(
	fortwrap
    fortwrapline
	freadahead
	fsubsplit
	fparamsplit
	fdeclsplit
	$fortran_extension
	$fortran_token
    $fortran_number
);
$VERSION = '0.01';

$fortran_token = '[_a-zA-Z][_\w]*';
$fortran_number =
	'[+-]?(?:\d+(?:.(?:\d+)?)?|.\d+)(?:[eEdD][+-]?\d+)?';
$fortran_extension = 'f';

1;

sub fortwrap {
    my ($textin) = @_;
    my $text = '';

    while ( $textin =~
           /([^\n]*\n     # Any stuff up to the teminating newline
             |            # or if there isn't a terminating newline
             .+)/gxs      # all other characters (but not the empty string)
           ) {
        my $line = $1;

        $line = fortwrapline( $line );

        $text .= $line;
    }

    return $text;
}

sub fortwrapline {
    my ($line) = @_;

    if ( length($line) > 72 ) {
        # This line could be too long, but isn't necessarily
        # if it includes the end-of-line character
        if ( length($line) > 73 || substr($line, 72, 72) ne "\n" ) {
            # The line is definitely too long
            my $rest = substr $line, 6;
            $line = substr $line, 0, 6;
            my $iscomment = $line =~ /^[*cC]/;
            my $term = chomp($rest) ? "\n" : '';

            while ( $rest ) {
                # Try to split at whitespace (or the end of string)
                $rest =~ s/^(.{0,66})(?:\s+|$)//
                ||  # Try to split at non-token characters
                $rest =~ s/^(.{0,65}[^_\w])//
                ||  # split any way you can
                $rest =~ s/^(.{66})//x
                || die "This can't happen";

                # At this point, rest is either empty, or contains
                # something other than whitespace.
                if ( $rest ) {
                    if ( $iscomment ) {
                        $line = $line . "$1\n*     ";
                    } else {
                        $line = $line . "$1\n     +";
                    }
                } else {
                    $line = $line . $1;
                }
            } # end while
            $line .= $term;
        } #end if the line is definitely too long
    } # end if the line may be too long

    return $line;
}

sub freadahead {
    my ($line, $FH, $filename) = @_;

    # First kill all leading/trailing whitespace and !-style comments
    for ($line) {
	s/!.*$//;
        s/^\s+//;
        s/\s+$//;
    }
    my $nextline = <$FH>;
    # Skip comment lines. (Does f77 allow comments in continuations?
    # Shug. My complier does, so I do too)
    while ($nextline && $nextline =~ /^[*cC]/) {
        $nextline = <$FH>;
    }
    while ( $nextline && $nextline =~ /^     \S/ ) {
        # This is a continuation line, append it to the
        # current line.
        # Pretty it up first, by removing leading/trailing
        # whitespace. It is *essential* that the trailing \n
        # be removed
        for ( $nextline ) {
            s/^     \S\s*/ /;
            s/\s+$//;
        }
        $line = $line . $nextline;
        $nextline = <$FH>;
        # Skip comment lines
        while ($nextline && $nextline =~ /^[*cC]/) {
            $nextline = <$FH>;
        }
    }
    return ($line, $nextline)
}

sub fparamsplit {
	my ($paramStr) = @_;
	my (@params, @ans);

	# Allow, but don't preserve, embedded whitespace
	$paramStr =~ s/\s+//g;

	$paramStr =~ s/^parameter\((.+)\)$/$1/
	or die "Syntax error in parameter declaration\n";

    @params = split /,/, $paramStr;
    @ans = ();
    foreach (@params) {
		/^($fortran_token)=(.+)$/o
			or die "Syntax error in parameter declaration\n";
		push @ans, $1, $2;
	}

    return @ans;
}

sub fsubsplit {
    my ($decl) = @_;

    # Allow, but don't preserve, embedded whitespace
    $decl =~ s/\s+//g;

    $decl =~
        /^subroutine
         ($fortran_token)  # The function name
         (?: \(            # The '(' starting the (optional) arg list
         ([^)]*)           # The arguments (anything that is not ')'
         \) )              # The close of the optional arg list
         $/xo;             # and nothing else

    my ($subname, $args) = ($1, $2);
    $subname or die "Syntax error in subroutine declaration\n";

    my @args = ();
    if ( $args ) {
        @args = split /,/, $args;
        foreach (@args) {
            /^$fortran_token$/
				or die "Syntax error in subroutine declaration\n";
        }
    }
    return ($subname, @args);
}

sub fdeclsplit {
	my ($line) = @_;

	my (@decl,@names);

	# Allow, but don't preserve, whitespace
	$line =~ s/\s//g;
	$line =~ s/\s//g;

	$line =~
		/^(integer |
		   doubleprecision |
		   character\*\d+ |
		   dimension)(.*)/ix
		or die "Line is not a fortran declaration\n";

	push @decl, $1, split /,/, $2;

	return @decl;
}

__END__
