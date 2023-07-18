#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Cwd;


# Converts abundance table to text file that can be converted to Krona graph with ktImportText
#
# input: abundance table
#	Table MUST be in tab-deleted format with a header row.
#	Table MUST include the following headers, any order is accepted, other headers are ignored:
#		taxid
#		readn
#		species
#		species_id
#		genus
#		genus_id
#		family
#		family_id
#		order
#		order_id 
#		class
#		class_id
#		phylum
#		phylum_id
#		superkingdom
#		superkingdom_id
#	--sample str	: name of the biological sample (optional)
#	--tail num	: number in percents [0,100], taxa with least abundance summing to this amount will be ignored
#	--rgtaxid	: taxid for the host genome. When given, this will be excluded from the tailed data.
#
# output: taxonomy profile in CAMI format (see https://github.com/bioboxes/rfc/blob/master/data-format/profiling.mkd)


my $usage=      "USAGE: $0 abund_table [-s|sample id -t|tail num --rgtaxid] 1> krona.txt\n".
		"\n".
		"Converts abundance table to Krona text format. Konvert further to Krona graph with ktImportText\n\n";
		
my $help 	= !1;
my $sampleid	= cwd();
my $tail	= 0;
my $rgtaxid	= 0;
GetOptions('help|h' 	=> \$help,
	'sample|s=s'	=> \$sampleid,
	'tail|t=f'	=> \$tail,
	'rgtaxid=i'	=> \$rgtaxid) or die $usage;
if((scalar @ARGV < 1) || $help){ die $usage; }


# ABUND TABLE <HEADER>
my @headers_taxnames = ('superkingdom','phylum','class','order','family','genus','species');
my @headers_taxids = ();
for my $h(@headers_taxnames){ push(@headers_taxids,join('_',$h,'id'))};
my @headers_required;
push(@headers_required,@headers_taxnames,@headers_taxids,'readn','taxid');
#my @headers_optional = ('strain');
my %headeri;
my $file = shift(@ARGV);
open(IN,"<$file") or die "Can\'t open $file: $!\n";
my $l=<IN>;
chomp($l);
$l =~ s/^#//;
my @headers= split(/\t/,$l,-1);
#DEBUG
#print "headers:$headers[0]\n";
for my $h(@headers_required){
	my $i = array_search(\@headers,$h);
	if($i < 0){ print STDERR "ERROR: $0: missing header: $h\n"; die $usage;}
	$headeri{$h} = $i;
}


# TAXA IN THE TAIL
my %tail_taxa;
if($tail > 0){
my %abund_species; 	# species_taxid > readn
my $readn_tot	= 0;
my $readn_host	= 0;
while($l=<IN>){
	if($l =~ m/^[#!]/){next;}
	chomp($l);
	my @sp		= split(/\t/,$l,-1);
	my $species_id  = $sp[$headeri{'species_id'}]; 
	my $abund	= $sp[$headeri{'readn'}];
	if($rgtaxid &&  ($rgtaxid eq $species_id)){
			$readn_host = $abund;
	}	
	$readn_tot	+= $sp[$headeri{'readn'}];
	if(!defined($abund_species{$species_id})){ $abund_species{$species_id} = 0;}
	$abund_species{$species_id} += $sp[$headeri{'readn'}];
}

my $cumsum = 0;
for my $k(sort {$abund_species{$a} <=> $abund_species{$b}} keys %abund_species){
	$cumsum += $abund_species{$k};
	if($cumsum > (($readn_tot-$readn_host)*$tail/100)){	# host reads excluded from tail cut-off
		last;
	}
	$tail_taxa{$k}= $abund_species{$k};
}
close(IN);
open(IN,"<$file") or die "Can\'t open $file: $!\n";
my $l=<IN>;
}



# print:READN TAB TAB-DELIMITED-TAXPATH


while($l=<IN>){
	if($l =~ m/^[#!]/){next;}
	chomp($l);
	my @sp = split(/\t/,$l,-1);
	# reorder data by header
	my %data;
	for my $h(@headers_required){
		$data{$h}	= $sp[$headeri{$h}];
	}
	if( $tail > 0){	# smalles taxa are ignored, e.g. the smallest 1%
		if(defined($tail_taxa{$data{'species_id'}})){
			next;
		}
	}	
	
	my @taxpathsn_line;
	for my $h(@headers_taxnames){
		if($data{$h} ne ''){
			push(@taxpathsn_line,$data{$h});
		}
	}
	
	# READN TAB-DELIMITED TAXPATH
	printf("%u\t%s\n", $data{'readn'}, join("\t",@taxpathsn_line));
}
close(IN);


sub array_search {
    my ($arr, $elem) = @_;
    my $idx= -1;
    for my $i (0..$#$arr) {
        if ($arr->[$i] eq $elem) {
            $idx = $i;
            last;
        }
    }
    return $idx;            
}






