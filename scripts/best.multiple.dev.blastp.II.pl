#!/usr/bin/perl
use strict;
use warnings;

## Goal
## Identify the best alignment(s) for each of the contigs 

  my $file = $ARGV[0] ;
  open(IN,"<$file") or die "Can\'t open $file: $!\n";
  my @INFO  =  <IN> ;

  my $output_file = $ARGV[1] ;   
  open(MP,">$output_file") or die "Can\'t open $output_file: $!\n";

   my $best_file = $ARGV[2] ;   
   open(BF,">$best_file") or die "Can\'t open $best_file: $!\n";
   

## args
## file : file output of blastn

  my %Mt ;
  my %Ct ;

# Field
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
     my $end_subject = $ar[9]	;
     my $line = join("\t",  @ar );
     
   # print  join("\t", "pre-filtration", $id, $target, $iden, $cov),"\n" ;    

     if(($iden >= 10) and ($cov >= 10)){

 # print  join("\t", "filtration", $id, $target, $iden, $cov),"\n" ;      
	 
      if($start_subject > $end_subject){
          
          my $tmp = $end_subject;
	  $end_subject = $start_subject;
          $start_subject = $tmp;
          
              }
   
## Proportion of the target sequence covered by the alignment

     my $proportion = (($end_subject-$start_subject)/$length_subject)*100;
     $Mt{$id}{$evalue}{$iden}{$cov}{$line} = "" ;
      

        }
}

  my @reads = (keys %Mt);

## For each of the contigs the alignments are ranked according to the evalue, identity, and cover of the query sequence in that order
 
while(my $rd  = shift @reads){

   my @nm_ord =  sort { $a  <=> $b  } (keys %{$Mt{$rd}});
   my $ev  = shift @nm_ord;

### Matrix ordered according to the e-value
### If the matrix has a single reference sequence this is the best hit  

   
   my @nm_por = sort { $b  <=> $a  }(keys %{$Mt{$rd}{$ev}});
   
   if ( scalar(@nm_por)  == 1){

     my $id  = shift @nm_por  ;
     my @Coverage = sort { $b  <=> $a  }(keys %{$Mt{$rd}{$ev}{$id}});
     my $cov =  shift @Coverage;
     my @Lines = (keys %{$Mt{$rd}{$ev}{$id}{$cov}});

##
  #   print  join("\t", "One e-value", $rd, $ev, $id, $cov),"\n" ;
##     
     my $line = shift  @Lines ;
     print BF $line,"\n";
     my @ar = split(/\t/, $line);
     my $ky = $ar[0] ;
     print MP join("\t", $ky, $ar[1] ),"\n" ; 
   
            }
### if there is more than one target sequence with the same e-value
### The ordered list is traversed from the highest identity to the lowest, within each of these values if there is more (or not) than one coverage related to the same identity value, it is evaluated if the difference is less than ten; Minor differences are stored in a table that is ordered from lowest to highest value to determine which coverage is closest to identity.
### In the event that there is no alignment in which the difference between identity and coverage is not less than ten, the one with the highest identity is taken (variables with the first prefix).
   
   else{

     my $ind = 0  ;
     my $first_cov = 0;
     my $first_id = 0;
     my $first_ind_id = 0;
     my $first_ind_cov = 0;
     
   while( ( my $id  = shift @nm_por ) and ( $ind != 1 )){

       
    my @Coverage = sort { $b  <=> $a  }(keys %{$Mt{$rd}{$ev}{$id}});


    my $size = scalar @Coverage;
    print  join("\t", "before", $rd, $ev, $id, $size),"\n" ;     



    ### hash where the key is the difference between the identity and the coverage and the value is the coverage

    if($first_ind_id == 0 ){
        $first_id = $id ;
	$first_ind_id = 1;
          }
    
      my %VL ; 
    
       while( my $cov =  shift @Coverage ){

#print  join("\t", $rd, $ev, $id, $cov),"\n" ; 
	   
        if($first_ind_cov == 0 ){

	  $first_cov = $cov  ;
	  $first_ind_cov = 1 ;

	        }

	   
	   my $diff = $id  -  $cov ;

	   if( $diff < 0){
               my $tmp =  $diff * (-1);
	       $diff = $tmp;
	           }

	   if($diff < 10){

	       $VL{$diff}  = $cov ;  
	       
	           }
                 }


    
    if(%VL){

	$ind = 1;
        my @cv_ord =  sort { $a  <=> $b  } (keys %VL);	
        my $dife = shift @cv_ord ;
	my $cov_best = $VL{$dife} ;
	my @Lines = (keys %{$Mt{$rd}{$ev}{$id}{$cov_best}});
	my $line = shift  @Lines ; 
	print BF $line,"\n";
	my @ar = split(/\t/, $line);
	my $ky = $ar[0] ;
	print MP join("\t", $ky, $ar[1] ),"\n" ; 

             }
           }

     if($ind == 0){

  # print  join("\t", $rd, $ev, $first_id, $first_cov),"\n" ; 
	 
     my @Lines = (keys %{$Mt{$rd}{$ev}{$first_id}{$first_cov}});
     my $line = shift  @Lines ; 
     print BF $line,"\n";
     my @ar = split(/\t/, $line);
     my $ky = $ar[0] ;
     print MP join("\t", $ky, $ar[1] ),"\n" ; 
      
      }
### else end
   }
}
  
          
close MP;
close BF;
