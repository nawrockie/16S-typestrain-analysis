use strict;

my $usage = "perl classify-given-seqstat.pl <map file> <seqstat file>";
if(scalar(@ARGV) != 2) { 
  die $usage;
}

my ($map_file, $seqstat_file) = (@ARGV);

my $nseq;
my %seqidx_H = ();
my %seqlen_H = ();
ribo_ParseSeqstatFile($seqstat_file, \$nseq, \%seqidx_H, \%seqlen_H);

open(MAP, $map_file) || die "ERROR unable to open $map_file for reading";
while(my $line = <MAP>) { 
  if($line !~ m/^\#/) { 
    if($line =~ /^N\S+\|(\S+)\:(\d+)\-(\d+)\|/) { 
      my $seqname = $1;
HERE HERE 

      my $length = $seqlen_H{


#################################################################
# Subroutine : ribo_ParseSeqstatFile()
# Incept:      EPN, Wed Dec 14 16:16:22 2016
#
# Purpose:     Parse an esl-seqstat -a output file.
#              
# Arguments: 
#   $seqstat_file:            file to parse
#   $nseq_R:                  REF to the number of sequences read, updated here
#   $seqidx_HR:               REF to hash of sequence indices to fill here
#   $seqlen_HR:               REF to hash of sequence lengths to fill here
#
# Returns:     Total number of nucleotides read (summed length of all sequences). 
#              Fills %{$seqidx_HR} and %{$seqlen_HR} and updates 
#              $$max_targetname_length_R, $$max_length_length_R, and $$nseq_R.
# 
# Dies:        If the sequence file has two sequences with identical names.
#              Error message will list all duplicates.
#              If no sequences were read.
#
################################################################# 
sub ribo_ParseSeqstatFile { 
  my $nargs_expected = 4;
  my $sub_name = "ribo_ParseSeqstatFile";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($seqstat_file, $nseq_R, $seqidx_HR, $seqlen_HR) = @_;

  open(IN, $seqstat_file) || die "ERROR unable to open esl-seqstat file $seqstat_file for reading";

  my $nread = 0;            # number of sequences read
  my $tot_length = 0;       # summed length of all sequences
  my $targetname;           # a target name
  my $length;               # length of a target
  my %seqdups_H = ();       # key is a sequence name that exists more than once in seq file, value is number of occurences
  my $at_least_one_dup = 0; # set to 1 if we find any duplicate sequence names

  while(my $line = <IN>) { 
    # = lcl|dna_BP331_0.3k:467     1232 
    # = lcl|dna_BP331_0.3k:10     1397 
    # = lcl|dna_BP331_0.3k:1052     1414 
    chomp $line;
    #print $line . "\n";
    if($line =~ /^\=\s+(\S+)\s+(\d+)/) { 
      $nread++;
      ($targetname, $length) = ($1, $2);
      if(exists($seqidx_HR->{$targetname})) { 
        if(exists($seqdups_H{$targetname})) { 
          $seqdups_H{$targetname}++;
        }
        else { 
          $seqdups_H{$targetname} = 2;
        }
        $at_least_one_dup = 1;
      }
        
      $seqidx_HR->{$targetname} = $nread;
      $seqlen_HR->{$targetname} = $length;
      $tot_length += $length;
    }
  }
  close(IN);
  if($nread == 0) { 
    die "ERROR did not read any sequence lengths in esl-seqstat file $seqstat_file, did you use -a option with esl-seqstat";
  }
  if($at_least_one_dup) { 
    my $i = 1;
    my $die_string = "\nERROR, not all sequences in input sequence file have a unique name. They must.\nList of sequences that occur more than once, with number of occurrences:\n";
    foreach $targetname (sort keys %seqdups_H) { 
      $die_string .= "\t($i) $targetname $seqdups_H{$targetname}\n";
      $i++;
    }
    $die_string .= "\n";
    die $die_string;
  }

  $$nseq_R = $nread;

  return $tot_length;
}
