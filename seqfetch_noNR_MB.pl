#! /usr/bin/perl

# usage: ./seqfetch_noNR_MB.pl output_file.fa
### This is to extract mixed random subseq of two organisms with user defined total number and ratio

use strict;
use warnings;

# User defined total number of mixed Mus/Bacterial sequences and their ratio
use constant mixedseq_count => 5000; 
use constant ratio_mb => 0.8; # Mus:Bacterial = 0.8


# get the respective Mus and Bacterial directories and output file name from command line 

my $mdirname = "/rsgrps/bhurwitz/hurwitzlab/data/reference/mouse_genome/20141111";
my $bdirname = "/rsgrps/bhurwitz/hurwitzlab/data/reference/mouse_genome/20141111/bacteria";


my $outfile = $ARGV[0];


my @mfiles = odir ($mdirname);
my @bfiles = odir ($bdirname);



# pass the subseq count, regular expression (file name), directory, file list of Mus and Bacteria to the subroutine seqfetch
# open all beginning with mm_alt_Mm_Celera_*.fa for the Mus directory
# open all beginning with SRR*.fasta for the Bacterial directory

my $mus_matcher = qr/^mm_alt_Mm_Celera_(\S+).fa$/; # such as mm_alt_Mm_Celera_chr1.fa and mm_alt_Mm_Celera_unplaced.fa
my $mus_dir = $mdirname;
my $mus_count = int(mixedseq_count * ratio_mb/22); #total 22 files in the Mus directory
my $mOUT = seqfetch($mus_count,$mus_matcher, $mus_dir, @mfiles);


my $bac_matcher = qr/^SRR(\d+).fasta$/; # such as SRR033547.fasta
my $bac_dir = $bdirname;
my $bac_count = (mixedseq_count - $mus_count)/3; #total 3 files in the Bacterial directory
my $bOUT = seqfetch($bac_count,$bac_matcher, $bac_dir, @bfiles);

print "Subseq count is $mus_count and $bac_count for individual Mus chr and for bacterial files, respectively\n";



##### subroutine odir: get a list of all files in a directory; ignore all files beginning with .                                                                  
sub odir {
    my ($dirname) = @_ ;
    opendir(my $dh, $dirname) || die "can't opendir $dirname: $!";
    my @files = grep { /^[^\.]/ && -f "$dirname/$_" } readdir($dh);
    closedir $dh;
    return @files;
}



##### subroutine seqfetch: extract qualified subseq

sub seqfetch {
  # global varibles and array 
  my $gi;
  my $seq;
  my $output;
  my $count;
  my @startarray;
  my @sorted_startarray;


  # loop through the files that match the defined pattern and write the output to Mus directory  
  my ($count, $match, $dirname, @files) = @_;
  foreach my $file (@files) {
    if ($file =~ /$match/) {
	 $gi = $1;
         open(DB,  "< $dirname/$file") or die "cant open $dirname/$file for reading";
         open(OUT, ">> $mdirname/outfile.fa") or die "cant open $mdirname/$1 for appending";


    my $output= 0;
    my $seqtot = 0;
    my $seq = "";
    while (<DB>) {   unless (/^>.+/) {
	chomp ($_);
        $seq .= $_; #single line for one seq header
                 }
    }
#print "header is $gi, seq is $seq\n"; #check if right header and seq retrieved. 

    while (1) {

        # if subseq count meets the requirement, then exit
	if ($output >= $count) {
          last;
        }

	my @sorted_startarray = ();
	my @startarray = ();
        my $seq_dif = $count - $output;
        while ($seq_dif >= 1) {
        # get random starting number
          my $seqtot = length($seq);
          my $startar = int(rand($seqtot));
          push(@startarray, $startar);  #append an element to the array, the element is a random number [0,9] if length is 10
          $seq_dif--;
          @sorted_startarray = sort {$a <=> $b} @startarray;
        }
#print "Start array is @sorted_startarray\n"; #check sorted random start array generated.
	
        # output the qualified subseq 
	foreach (@sorted_startarray) {
	  chomp ($_);
          my $start = $_ - 1;
          my $length = 200; #length of subseq
	  my $end   = $start + $length;
    	  my $outseq= substr ($seq, $start, $length);
          # disqualify the subseq with length < 200, or including 'N'
	      next if length($outseq) < $length;
                  next if $outseq=~ /N/;  
                      print OUT ">$gi\_$start\_$end\n", "$outseq\n";
		      ++$output;
        }   
    }

print "$output sequences have been exported for $gi\n";

      close(DB); 
    } else {
      print "$file didn't match\n";
      }
  }
}
