#! /usr/bin/perl -w -s

my @seqs;

my $header="";
my $h1="";
my $seq="";
my $output;
#my $cutoff=75;

while (<STDIN>) {
	chomp;
	if ($_=~/\>/) {
		$h1=$_;
		if ($seq) {
			$seq=dna_to_protein($seq);
			push @seqs,([$header, $seq]);
		}
		#push @seqs,([length($seq),$header, $seq ]);
		$seq="";		
	} else{
		$header=$h1;
		$seq.=$_;
	}
}

$seq=dna_to_protein($seq);
if($header) {
	push @seqs,([$header,$seq]) if ($seq ne "");# && (length($seq)>=$cutoff);
} else {
	push @seqs,([">unknown",$seq])
}	
	
#push @seqs,([length($seq),$header, $seq ]);

#shift @seqs;

my $counter=1;
foreach(@seqs) {
	$output.= sprintf "%s\n%s\n",@{$_};
}
		
syswrite STDOUT, $output;

sub dna_to_protein {
#loops through a nucleotide sequence and feeds three character codons to AminoOut function, returns protein chain.



	my ($dna) = @_;
	my $protein = "";

	my %codon2aa = qw(
	    TCA  S  TCC  S  TCG  S  TCT  S  TTC  F  TTT  F  TTA  L  TTG  L
	    TAC  Y  TAT  Y  TAA  *  TAG  *  TGC  C  TGT  C  TGA  *  TGG  W
	    CTA  L  CTC  L  CTG  L  CTT  L  CCA  P  CCC  P  CCG  P  CCT  P
	    CAC  H  CAT  H  CAA  Q  CAG  Q  CGA  R  CGC  R  CGG  R  CGT  R
	    ATA  I  ATC  I  ATT  I  ATG  M  ACA  T  ACC  T  ACG  T  ACT  T
	    AAC  N  AAT  N  AAA  K  AAG  K  AGC  S  AGT  S  AGA  R  AGG  R
	    GTA  V  GTC  V  GTG  V  GTT  V  GCA  A  GCC  A  GCG  A  GCT  A
	    GAC  D  GAT  D  GAA  E  GAG  E  GGA  G  GGC  G  GGG  G  GGT  G
	);

	while ( $dna =~ /(...)/g ) {
	    exists $codon2aa{ $1 } or $codon2aa{ $1 }="";#die qq[Bad codon "$1"!!\n];
	    $protein .= $codon2aa{ $1 };
	}

	return $protein;
}	


sub AminoOut {
#converts a three character nucleotide sequence to it's specified Amino Acid
#
#NOTE - hash copied verbatim from one of the practicals.
		
	my ($codon) =@_;
	$codon = uc($codon);
	my (%codon_dictionary) = (
		GCT => 'A', GCC => 'A', GCA => 'A', GCG => 'A',  
		TGT => 'C', TGC => 'C',
		GAT => 'D', GAC => 'D',
		GAA => 'E', GAG => 'E',
		TTT => 'F', TTC => 'F',
		GGT => 'G', GGC => 'G', GGA => 'G', GGG => 'G',
		CAT => 'H', CAC => 'H',
		ATT => 'I', ATC => 'I', ATA => 'I',
		AAA => 'K', AAG => 'K', 
		TTA => 'L', TTG => 'L', CTT => 'L', CTC => 'L', CTA => 'L', CTG => 'L',
		ATG => 'M', 
		AAT => 'N', AAC => 'N',
		CCT => 'P', CCC => 'P', CCA => 'P', CCG => 'P', 
		CAA => 'Q', CAG => 'Q', 
		CGT => 'R', CGC => 'R', CGA => 'R', CGG => 'R', AGA => 'R', AGG => 'R',
		TCT => 'S', TCC => 'S', TCA => 'S', TCG => 'S', AGT => 'S', AGC => 'S',
		ACT => 'T', ACC => 'T', ACA => 'T', ACG => 'T',
		GTT => 'V', GTC => 'V', GTA => 'V', GTG => 'V', 
		TGG => 'W', TGA => '*',
		TAT => 'Y', TAC => 'Y',
		TAA => '*', TAG => '*',
	);

	if (exists ($codon_dictionary{$codon})) {
		return $codon_dictionary{$codon};
	} else {
		#print STDERR "Bad codon: ";
		#print STDERR "$codon!!\n\n";
		return "-";
	}
}


