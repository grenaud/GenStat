#!/usr/bin/perl


use strict;
use warnings;
use Data::Dumper;
use Cwd 'abs_path';
use File::Basename;





my $pathtoperlsc = dirname(abs_path($0))."/";




sub formatDivergence{
  my ($toprint)=@_;
  return sprintf("%.3f", (100*$toprint))."";
}


sub computeBoot{
  my ($shortBranch,$commonBranch) = @_;
  my $cmd=$pathtoperlsc."binom.R $shortBranch $commonBranch";
  #die $cmd;
  my @output=split("\n",`$cmd`);
  chomp($output[0]);
  chomp($output[1]);

  return (formatDivergence($output[0]),
	  formatDivergence($output[1]));
}

#die Dumper(computeBoot(19,91));

if($#ARGV == -1){
  die "Usage ./combineOutputs.pl [output1] [output2] ....\n"."this script merely combines the output from many divergence outputs into one\n";

}


my $columnsPerCategory=10;
my $dummy=-20;
my $columns=$dummy;
my $arraySum=[];

for my $arg (@ARGV){
#  print "arg ".$arg."\n";
  open(INFILE,$arg) or die "Cannot open file $arg";
  while(my $line=<INFILE>){
    chomp($line);

    if($line=~/^\s*$/ ||
       $line=~/^#/  ){
      next;
    }

    my @arraytemp=split("\t",$line);

    if($arraytemp[1] eq "N/A"){
      next;
    }

    my $numberOfColumns = $#arraytemp;
    #print "numberOfColumns $numberOfColumns\n";
    if( ($numberOfColumns % $columnsPerCategory ) != 0 ){
      die "Line $line in file ".$arg." should have a multiple of $columnsPerCategory as the number of columns, found $numberOfColumns\n";
    }

    if($columns ==  $dummy){
      $columns = $numberOfColumns;
      for(my $i=0;$i<($numberOfColumns /$columnsPerCategory );$i++){
	push(@{$arraySum},[0,0,0,0]);
      }
    }else{
      if($columns !=  $numberOfColumns){
	die "Line $line should have the same number of columns as previous files (".$columns."), found $numberOfColumns\n";
      }
    }


    for(my $i=0;$i<($numberOfColumns /$columnsPerCategory );$i++){
      my $indexBase=$i*$columnsPerCategory+1;
      for(my $j=0;$j<4;$j++){
	@{@{$arraySum}[$i]}[$j]+=$arraytemp[$indexBase+$j];
      }
    }
  }
  close(INFILE);
}

#0 no mutation
#1 common
#2 reference
#3 sample

my $lastIndex=($columns /$columnsPerCategory );
for(my $i=0;$i<$lastIndex;$i++){
  for(my $j=0;$j<4;$j++){
    my $valprint=@{@{$arraySum}[$i]}[$j];
    print $valprint."\t";
  }

  #R/(R+C)
  if( (@{@{$arraySum}[$i]}[2]+@{@{$arraySum}[$i]}[1]) != 0){
    print join("\t",(formatDivergence(@{@{$arraySum}[$i]}[2]/(@{@{$arraySum}[$i]}[2]+@{@{$arraySum}[$i]}[1])),
		     computeBoot( @{@{$arraySum}[$i]}[2], @{@{$arraySum}[$i]}[1] ) ) );
		     #formatDivergence(lowerConf(@{@{$arraySum}[$i]}[2],@{@{$arraySum}[$i]}[1])),
		     #formatDivergence(highConf( @{@{$arraySum}[$i]}[2],@{@{$arraySum}[$i]}[1]))) );
  }else{
    #print "NaN\tNan\tNan\t";
    print join("\t",("NaN","NaN","NaN"));
  }
  #if($i < ($lastIndex-1)){
  print "\t";
  #}

  #S/(S+C)
  if(   (@{@{$arraySum}[$i]}[3]+@{@{$arraySum}[$i]}[1])  ){
    print join("\t",(formatDivergence(@{@{$arraySum}[$i]}[3]/(@{@{$arraySum}[$i]}[3]+@{@{$arraySum}[$i]}[1])),
		     computeBoot( @{@{$arraySum}[$i]}[3], @{@{$arraySum}[$i]}[1] ) ) );

		     #formatDivergence(lowerConf(@{@{$arraySum}[$i]}[3],@{@{$arraySum}[$i]}[1])),
		     #formatDivergence(highConf(@{@{$arraySum}[$i]}[3],@{@{$arraySum}[$i]}[1]))) );
  }else{
    #print "NaN\tNan\tNan\t";
    print join("\t",("NaN","NaN","NaN"));
  }

  if (  $i != ($lastIndex-1)){
    print "\t";
  }

}
print "\n";
