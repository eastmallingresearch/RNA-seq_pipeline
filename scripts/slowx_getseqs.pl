#! /usr/bin/perl -w -s
use Time::HiRes qw( time );

# usearch bug in fastx_getseqs - this is my replacement for this program...
# Requires a lot of memory ~ 32G

my $file = $ARGV[0];
my $output;

my %headers;
my $h1="";
my $seq="";

my $start = time();

print STDERR "Reading Headers...\n";

while (<STDIN>) {
	chomp;
	$headers{$_}="";	
}
my $end = time();

my $pr = sprintf("%.2f\n", $end - $start);
$start = $end;

print STDERR "Headers read in $pr\n";
print STDERR "Reading fastq...\n";

open(FASTX,$file);

#my $skip

while (<FASTX>) {
	chomp;
#	print STDERR "$.\n";
	if ($.%4==1) {
		my $header=substr $_, 1;
		if (exists $headers{$header}) {
			$output.="$_\n";
			$output.=<FASTX>;
			$output.=<FASTX>;
			$output.=<FASTX>;
			$headers{$header}=$output;
			#syswrite STDOUT, $output;
			$output="";
		}
	}
}

$end = time();
$pr = sprintf("%.2f\n", $end - $start);

print STDERR "fastq processed in $pr\n";

print STDERR "writing output...\n";
$start = $end;

foreach my $key (sort keys %headers)  {
	print"$headers{$key}";
}

$end = time();
$pr = sprintf("%.2f\n", $end - $start);
print STDERR "Written in $pr\n";


#syswrite STDOUT, $output;
