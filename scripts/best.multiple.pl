#!/usr/bin/perl
use strict;
use warnings;

## Goal
## Identify the best alignment(s) for each of the contigs 

  my $file = $ARGV[0] ;
  open(IN,"<$file") or die "Can\'t open $file: $!\n";
  my @INFO  =  <IN> ;

## args
## file : file output of blastn

  my %Mt ;
my %Ct ;

# Fields
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send  evalue  bitscore qcovs slen sgi  staxids sscinames stitl


  while( my $ln =  shift  @INFO ){

     chomp($ln);
     my @ar = split(/\t/, $ln);
     my $id = $ar[0] ;        
     my $target = $ar[1] ;  
     my $iden = $ar[2] ;
     my $cov = $ar[12] ;
     my $length_subject= $ar[13] ;
     my $evalue = $ar[10] ; 
     my $start_subject = $ar[8] ;
     my	$end_subject = $ar[9]	;
  
      if($start_subject > $end_subject){
          
          my $tmp = $end_subject;
	  $end_subject = $start_subject;
          $start_subject = $tmp;
          
              }
   
## Proportion of the target sequence covered by the alignment

     my $proportion = (($end_subject-$start_subject)/$length_subject)*100;
     my $new = join("\t", $ln, $proportion );  
	 
     $Mt{$id}{$evalue}{$iden}{$cov}{$new} = "" ;
      

     }
     

  my @reads = (keys %Mt);

## For each of the contigs the alignments are ranked according to the evalue, identity, and cover of the query sequence in that order
 
while(my $rd  = shift @reads){

   my @nm_ord =  sort { $a  <=> $b  } (keys %{$Mt{$rd}});
   my $id  = shift @nm_ord;
   my @nm_por = sort { $b  <=> $a  }(keys %{$Mt{$rd}{$id}});
   my $cv  = shift @nm_por  ;
   my @Coverage = sort { $b  <=> $a  }(keys %{$Mt{$rd}{$id}{$cv}});
   my $cov =  shift @Coverage;
   my @Lines = (keys %{$Mt{$rd}{$id}{$cv}{$cov}});

   my %nw  ;
   my $ky = "NA" ;


## Identified the best or all the best targets  
   
     while(my $line = shift @Lines){

	 my @ar = split(/\t/, $line);
	 $nw{$ar[1]} = "";
	 $ky = $ar[0] ;

            }


   my @info = keys %nw ;
   my $imp = join(",", @info); 

   print join("\t", $ky, $imp ),"\n" ;


     }
          
