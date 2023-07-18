#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Adds taxonomy information to the end of blast (or other) search result table
#

my $usage=      "USAGE: $0 -h -v blast.out nodes.dmp names.dmp species[,genus,family,order]\n".
		"\n".
		"blast.out   \t: blast result table, last col is assumed to be taxid\n".
		"            \t: if the first line start with \#, assumes this is a header line\n".
                "nodes.dmp   \t: node file from NCBI taxonomy\n".
                "names.dmp   \t: name file from NCBI taxonomy\n".
		"species[..] \t: comma-sep list of taxinfo to be printed:\n".
		"\t: species \t: species name\n".
		"\t: genus   \t: genus name\n".
		"\t: family  \t: family name\n".
		"\t: order   \t: order name\n".
		"\t: species_id\t: species tax id\n".
		"\t: etc.\n".
		"-i|ids      \t: print taxids after each taxon name\n".
		"-v           \t: print warning messages\n".
		"-h           \t: blast.out has headers\n".
                "\n".
		"nodes.dmp & names.dmp can be downloaded from:\n".
		"ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz\n";
my $headers 	= !1;
my $verb 	= !1;
my $printids	= !1;
GetOptions('headers|h' => \$headers,
		'verbal|v' => \$verb,
		'ids|i'	 => \$printids) or die $usage;

if(scalar @ARGV<4) { die $usage; }

my $blast_file		= shift(@ARGV);
my $node_dmp_file	= shift(@ARGV);
my $name_dmp_file	= shift(@ARGV);
my $format		= shift(@ARGV);
my @format_ranks	= split(/,/,lc($format));


# read taxonomy nodes
print STDERR "# reading $node_dmp_file\n";
my %taxnodes;
open(IN,"<$node_dmp_file") or die "Can\'t open $node_dmp_file: $!\n";
my $ln=0;
while(my $l=<IN>){
	$ln++;
	chomp($l);
	my @sp= split(/\t\|\t/,$l);
	if(scalar(@sp)<3){
		die "ERROR: missing data at $ln:$node_dmp_file\n";
	}
	$taxnodes{$sp[0]}= join("\t",$sp[0],$sp[1],$sp[2]);
}

# read taxonomy node names
print STDERR "# reading $name_dmp_file\n";
my %taxnames;
open(IN,"<$name_dmp_file") or die "Can\'t open $name_dmp_file: $!\n";
$ln=0;
while(my $l=<IN>){
	$ln++;
	chomp($l);
	my @sp= split(/\t\|\t/,$l);
	for(my $j=0; $j<scalar(@sp); $j++){
		$sp[$j] =~ s/^(\|\t)|(\t\|)$//g; # removing leading and trailing \t|\t characters
	}
	if(scalar(@sp)<4){
		die "ERROR: missing data at $ln:$name_dmp_file\n";
	}
	if( !defined($taxnames{$sp[0]}) ){ # if not defined we add any name available
		$taxnames{$sp[0]}= $sp[1];
	}
	elsif( lc($sp[3]) eq 'scientific name' ){ # else we replace saved name with scientific
		$taxnames{$sp[0]}= $sp[1];
	}
}

# USING %taxnodes and %taxnames hashes:
#
#my $id= 1656;
#my @sp= split(/\t/, $taxnodes{$id});
#print "hashnode: ",$taxnodes{$id},"\n";
#print "# node $id:\n";
#print "id         :",$sp[0],"\n";
#print "parent id  :",$sp[1],"\n";
#print "rank       :",$sp[2],"\n";
#print "name       :",$taxnames{$sp[0]},"\n";

# USING get_parent_taxids():
#
#print "# parents($id):\n";
#my @parent_taxids= get_parent_taxids(\%taxnodes,$id);
#foreach my $i(@parent_taxids){
#	my @sp = split(/\t/,$taxnodes{$i});
#	print "id   \t:$sp[0]\n";
#	print "rank \t:$sp[2]\n";
#	print "name \t:",$taxnames{$sp[0]},"\n";
#}
#exit(1);

print STDERR "# reading taxids from the last column in $blast_file\n";
open(IN,"<$blast_file") or die "Can\'t open $blast_file: $!\n";
$ln	= 0;
my $taxids_linked 	= 0;
my $taxids_empty	= 0;
my $taxids_unknown	= 0;
while(my $l=<IN>){
	$ln++;
	chomp($l);
	if($l =~ /^\s*$/){next;}
	
	# print headers if present, add headers for added taxonomy cols
	if( ($ln==1) && $headers ){ 
		print "$l";
		foreach my $taxon(@format_ranks){
			print "\t$taxon";
			if($printids){
			print "\t$taxon","_id";
			}
		}
		print "\n";
		next;
	}
	
	print "$l";
	my @sp		= split(/\t/,$l,-1);
	my @sp2		= split(/;/,$sp[$#sp]);
	my $hit_taxid	= '';
	if(scalar(@sp2)>0){$hit_taxid = shift(@sp2);};
	if( $hit_taxid eq ''){	# this is normal for SANS/BLAST tables with empty annotations (= no hit)
		for(my $j=0; $j<scalar(@format_ranks); $j++){
			print "\t";
			if($printids){
				print "\t";
			}
		}
		print "\n";	
		$taxids_empty++;
		# DEBUG print STDERR "empty taxid: $l\n";
		next;
	}
	elsif(!defined($taxnodes{$hit_taxid})){
		# skipping line
		if($verb){
			print STDERR "ERROR\tskipping unknown taxid\t$hit_taxid\n";
		}
		for(my $j=0; $j<scalar(@format_ranks); $j++){
			print "\t";
			if($printids){
				print "\t";
			}
		}
		print "\n";
		$taxids_unknown++;
		next;
	}
	
	
	@sp		= split(/\t/,$taxnodes{$hit_taxid},-1);
	#if( $sp[2] !~ m/species/i ){
	#	if($verb){
	#		print STDERR "WARNING: unexpected taxrank=$sp[2] for taxid=$sp[0] in $ln:$blast_file\n";
	#	}
	#}
	
	my @parent_taxids= get_parent_taxids(\%taxnodes,$hit_taxid);
	foreach my $format_rank(@format_ranks){
		my $rank_found	= 0;
		foreach my $p(@parent_taxids){
			@sp = split(/\t/,$taxnodes{$p},-1);
			if($sp[2] eq $format_rank){
				print "\t",$taxnames{$p};
				if($printids){
				print "\t",$p;
				}
				$rank_found = 1;
				last;
			}
		}
		if(!$rank_found){
			print "\t";
			if($printids){
				print "\t";
			}
		}
	}
	$taxids_linked++;
	print "\n";
}
my $taxids_processed	= ($taxids_empty + $taxids_unknown + $taxids_linked);
printf STDERR ("# taxids linked:\t%u/%u (%.2f%%)\n", $taxids_linked,$taxids_processed, ($taxids_linked/$taxids_processed*100));
if($taxids_empty > 0){
printf STDERR ("# taxids empty: \t%u/%u (%.2f%%)\n", $taxids_empty,$taxids_processed, ($taxids_empty/$taxids_processed*100));}
if($taxids_unknown > 0){
printf STDERR ("# taxids unknown:\t%u/%u (%.2f%%)\n", $taxids_unknown,$taxids_processed, ($taxids_unknown/$taxids_processed*100));}

# Returns parent tax ids, including the node itself		
sub get_parent_taxids{
	my $taxnode_ref	= shift(@_);
	my $taxid	= shift(@_);
	my @parent_taxids;
	push(@parent_taxids,$taxid);
	
	while( defined($taxnode_ref->{$taxid}) ){
		my $node	= $taxnode_ref->{$taxid};
		my @sp		= split(/\t/,$node,-1);
		if($sp[0] == $sp[1]){ 		# node is it's own parent
			last;
		}
		$taxid		= $sp[1];	# parent taxid
		push(@parent_taxids,$taxid);
		if($taxid == 1){		# root node
			last;
		}
	}
	return @parent_taxids;
}

