#!/usr/bin/perl


use strict;
use warnings;
use Data::Dumper;


sub formatDivergence{
  my ($toprint)=@_;
  return sprintf("%.3f", (100*$toprint))."%";
}


sub lowerConf {
  my ($shortBranch,$commonBranch) = @_;

  my $n=$shortBranch+$commonBranch;
  my $nsqrt = sqrt ( $n );
  #http://en.wikipedia.org/wiki/Checking_whether_a_coin_is_fair
  my $z = 1.9599;
  return ($shortBranch/$n)-($z/(2*$nsqrt));
}

sub highConf {
  my ($shortBranch,$commonBranch) = @_;

  my $n=$shortBranch+$commonBranch;
  my $nsqrt = sqrt ( $n );
  #http://en.wikipedia.org/wiki/Checking_whether_a_coin_is_fair
  my $z = 1.9599;
  return ($shortBranch/$n)+($z/(2*$nsqrt));
}


if($#ARGV == -1){
  die "Usage ./combineOutputs.pl [output1] [output2] ....\n"."this script merely combines the output from many divergence outputs into one\n";

}


my $columnsPerCategory=10;
my $dummy=-20;
my $columns=$dummy;
my $arraySum=[];

for my $arg (@ARGV){
  open(INFILE,$arg) or die "Cannot open file $arg";
  while(my $line=<INFILE>){
    chomp($line);
    my @arraytemp=split("\t",$line);
    if($arraytemp[1] eq "N/A"){
      next;
    }

    my $numberOfColumns = $#arraytemp;

    if( ($numberOfColumns % $columnsPerCategory ) != 0 ){
      #warn $line."\n";
      #for(my $col=0;$col<=($numberOfColumns);$col++){
      #warn $col."\t".$arraytemp[$col]."\n";
      #}
      die "Line $line in file ".$arg." should have a multiple of $columnsPerCategory as the number of columns, found $numberOfColumns\n";
    }

    if($columns ==  $dummy){
      $columns = $numberOfColumns;
      for(my $i=0;$i<($numberOfColumns /$columnsPerCategory );$i++){
	#print "i=$i\n";
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
		     formatDivergence(lowerConf(@{@{$arraySum}[$i]}[2],@{@{$arraySum}[$i]}[1])),
		     formatDivergence(highConf( @{@{$arraySum}[$i]}[2],@{@{$arraySum}[$i]}[1]))) );
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
		     formatDivergence(lowerConf(@{@{$arraySum}[$i]}[3],@{@{$arraySum}[$i]}[1])),
		     formatDivergence(highConf(@{@{$arraySum}[$i]}[3],@{@{$arraySum}[$i]}[1]))) );
  }else{
    #print "NaN\tNan\tNan\t";
    print join("\t",("NaN","NaN","NaN"));
  }

  if (  $i != ($lastIndex-1)){
    print "\t";
  }

}
print "\n";
