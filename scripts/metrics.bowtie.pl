#!/usr/bin/perl
use strict;
use warnings;

## Goal                                                                                               
## Get the number of paired-end and unpaired reads aligned to the sequences  


     my $file = $ARGV[0] ;
     open(IN,"<$file") or die "Can\'t open $file: $!\n";
     my @INFO  =  <IN> ; 

## Arguments
## file: bowtie2 output report file  

    my %Mt ;
    my $pre = "";

   while( my $ln  =  shift @INFO ){

         chomp($ln);

## Number total single reads

	if($ln =~  /.*reads; of these:.*/){

              $ln =~ s/\s*([0-9]+)\s+reads; of these:.*/$1/; 
              $Mt{total}{total} = $ln;                     

	      }
       elsif($ln =~  /.+were paired; of these:/){ 
             $pre = "paired" ;  
	        }

##  Number of paired-end reads (need to multiply by 2)

       elsif($ln =~ /.+aligned concordantly exactly 1 time.*/){

           $ln =~ s/\s+([0-9]+)\s+\([0-9.]+%\) aligned concordantly exactly 1 time.*/$1/;
           my $subf = "concordantly_1" ;
           $Mt{$pre}{$subf} = $ln;

	      }

##  Number of paired-end reads (need to multiply by 2) 

       elsif($ln =~ /.+aligned concordantly >1 times.*/){  

           $ln =~ s/\s+([0-9]+)\s+\([0-9.]+%\) aligned concordantly >1 times.*/$1/;
	   my $subf = "concordantly_>_1" ;
           $Mt{$pre}{$subf} = $ln;

             }

       elsif($ln =~  /.+pairs aligned concordantly 0 times; of these:.*/){
    	   $pre = "paired" ;
	        }

##  Number of paired-end reads (need to multiply by 2) 

       elsif($ln =~ /.*aligned discordantly 1 time.*/){

	  $ln =~ s/\s+([0-9]+)\s+\([0-9.]+%\) aligned discordantly 1 time.*/$1/;
          my $subf = "discordantly_1" ; 
	  $Mt{$pre}{$subf} = $ln;
  
                  }

##  Number of paired-end reads (need to multiply by 2) 
       
     elsif($ln =~ /.*aligned discordantly >1 times.*/){

          $ln =~ s/\s+([0-9]+)\s+\([0-9.]+%\) aligned discordantly >1 times.*/$1/;
          my $subf = "discordantly_>1" ;
          $Mt{$pre}{$subf} = $ln;

	 }

## only one of the reads from the paired-end line

    elsif($ln =~ /.*pairs aligned 0 times concordantly or discordantly; of these:*/){
         $pre = "paired_only_a_read_aligned" ;
  	    }

## Number single reads   

     elsif($ln =~ /.*aligned exactly 1 time.*/){

           $ln =~ s/\s+([0-9]+)\s+\([0-9.]+%\) aligned exactly 1 time.*/$1/;
	   my $subf = "aligned_1" ;
	   $Mt{$pre}{$subf} = $ln;
 
         }

## Number single reads 

     elsif($ln =~ /.*aligned >1 times.*/){

	  $ln =~ s/\s+([0-9]+)\s+\([0-9.]+%\) aligned >1 times/$1/;
          my $subf = "aligned_>1" ; 
	  $Mt{$pre}{$subf} = $ln;

           }


     elsif($ln =~ /.*were unpaired; of these:.*/){
         $pre = "unpaired" ;
             }
 
## Number single reads  

     elsif($ln =~ /.*aligned exactly 1 time/){

	 $ln =~ s/\s+([0-9]+)\s+\([0-9.]+%\) aligned exactly 1 time.*/$1/;
         my $subf = "exactly_1" ;
         $Mt{$pre}{$subf} = $ln;         

	 }

## Number single reads 

    elsif($ln =~ /.*aligned >1 times.*/){

	$ln =~ s/\s+([0-9]+)\s+\([0-9.]+%\) aligned >1 times.*/$1/;
        my $subf = "aligned_>1" ;
        $Mt{$pre}{$subf} = $ln;

           }

     }


  my @Type = keys %Mt;

  while(my $tp = shift @Type){
      my @cate = keys %{$Mt{$tp}};
     while(my $ct = shift @cate){ 
	 my $nm = $Mt{$tp}{$ct} ;
	 print join("\t",  $tp, $ct, $nm  ),"\n" ;
           }
       }

