#! /usr/bin/perl -w -s
use List::MoreUtils 'first_index';
use List::MoreUtils qw(indexes);


###########################################################
#
# Outputs all ORFs (longer than min_length) found in input sequence
# Optionally converts to protein and sorts by length (not implemented)
#
###########################################################

my @seqs;
my $min_length=75; #$ARGV[0] // 75;
#my $translate=$ARGV[1] // 0;
#my $sortme=$ARGV[2] // 0;

my $header="";
my $seq="";
my $h1="";
my @idx;
my $output;
my $cds;


while (<STDIN>) {
	chomp;
	if ($_=~/\>/) {
		$h1=(split / /,$_)[0];
		get_cds($seq,$h1) if $seq;
		$seq="";
	} else{
		$_=~s/\s+$//;
		$header=$h1;
		$seq.=$_;
	}
}

get_cds($seq,$h1);

foreach(@seqs) {
	$output.= sprintf "%s_%s_%s\n%s\n",@{$_};
}

syswrite STDOUT, $output;

sub get_cds {
	my ($s,$header) =  @_;
	my $rcs = reverse_compliment($s);
	my @ss;
	get_ss($s,$header,"1_F");
	get_ss($rcs,$header,"1_RC");
	$s =~ s/^.//s;
	$rcs =~ s/^.//s;
	get_ss($s,$header,"2_F");
	get_ss($rcs,$header,"2_RC");
	$s =~ s/^.//s;
	$rcs =~ s/^.//s;
	get_ss($s,$header,"3_F");
	get_ss($rcs,$header,"3_RC");
}

sub reverse_compliment {
	my ($s)=@_;
	$s=~tr/atcgATCG/tagcTAGC/;
	my $s2=reverse $s;
	return($s2);
}

sub get_ss {
	my ($s,$h,$class)=@_;
	my @ind = indexes { /TAG|TAA|TGA/ } ( $s =~ m/.../g ); # this finds all stop codons - which is also a map of potential ORFs
#foreach(@ind) {
#	print "$_\n";
#}
#exit;
	my $c=1;
	if (!@ind) {
		push @seqs, ([$h,$class,$c,$s]);
	} else {
		unshift(@ind,0) if $ind[0]!=0;
		push(@ind,length($s));
		while(scalar(@ind)>1) {
			my $l = substr($s,($ind[0]*3),($ind[1]*3 - $ind[0]*3));
			push @seqs, ([$h,$class,$c,$l]) if length($l)>=$min_length;
			shift(@ind);
			shift(@ind);
			$c++;
		}
			
	}
}

