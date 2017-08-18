#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

#require "epn-options.pm";

my $usage = "perl check-map-with-cmsearch.pl <map file> <16S-typestrain sequence file> <source sequence file> <file with list of known introns> <CM file> <output root>";
if(scalar(@ARGV) != 6) { 
  die $usage;
}

my ($map_file, $type_file, $source_file, $intron_file, $cm_file, $out_root) = (@ARGV);

# run esl-seqstat on both sequence files
my $type_seqstat   = remove_dir_path($type_file) . ".seqstat";
my $source_seqstat = remove_dir_path($source_file) . ".seqstat";

#run_command("esl-seqstat -a $type_file   > $type_seqstat",   0);
#run_command("esl-seqstat -a $source_file > $source_seqstat", 0);

my $type_nseq;
my %ltype_seqidx_H  = (); # key is long type name, value is sequence index in short sequence file
my %ltype_len_H  = (); # key is long type name, value is sequence length in short sequence file
parse_seqstat_file($type_seqstat, \$type_nseq, \%ltype_seqidx_H, \%ltype_len_H);

my %stype_seqidx_H = ();
my %stype_len_H = ();
my %stype2ltype_H   = ();
my %ltype2stype_H   = ();
shorten_keys(\%ltype_seqidx_H, \%ltype_len_H, \%stype_seqidx_H, \%stype_len_H, \%stype2ltype_H, \%ltype2stype_H);

my $source_nseq;
my %lsource_seqidx_H = ();
my %lsource_len_H = ();
parse_seqstat_file($source_seqstat, \$source_nseq, \%lsource_seqidx_H, \%lsource_len_H);

my %ssource_seqidx_H = ();
my %ssource_len_H = ();
my %ssource2lsource_H   = ();
my %lsource2ssource_H   = ();
shorten_keys(\%lsource_seqidx_H, \%lsource_len_H, \%ssource_seqidx_H, \%ssource_len_H, \%ssource2lsource_H, \%lsource2ssource_H);

my $line;
my $stype_name;   # short type name
my $ssource_name; # short source name
my $mapstart;
my $mapstop;
my $strand;
my $ltype_name;
my $lsource_name;
my %stype_ssource_H   = (); # key is short type name, value is short source name
my %ssource_stype_HA  = (); # key is short source name, value is array of type names that map to this source
my %stype_mapstart_H = (); # key is short type name, value is start index from map file
my %stype_mapstop_H  = (); # key is short type name, value is stop index from map file
my %stype_strand_H   = (); # key is short type name, value is strand from map file ('+' or '-')

open(MAP, $map_file) || die "ERROR unable to open $map_file for reading";
my $map_idx = 1;
my %map_stype_idx_H = (); # key is stype name, value is index in map file
my %map_stype_len_H = (); # key is stype name, value is length in map file (stop-start+1)
while(my $line = <MAP>) { 
  chomp $line;
  if($line !~ m/^\#/) { 
    if($line =~ /^(N\S+)\|(\S+)\.\d+\:(\d+)\-(\d+)\|/) { 
      ($stype_name, $ssource_name, $mapstart, $mapstop) = ($1, $2, $3, $4);
    }
    elsif($line =~ /^(N\S+)\|(\S+)\:(\d+)\-(\d+)\|/) { 
      ($stype_name, $ssource_name, $mapstart, $mapstop) = ($1, $2, $3, $4);
    }
    else { 
      die "ERROR unable to parse map line: $line";
    }
    if(! exists $stype2ltype_H{$stype_name}) { 
      die "ERROR type name $stype_name does not have a long name"; 
    }
    if(! exists $ssource2lsource_H{$ssource_name}) { 
      die "ERROR source name $ssource_name does not have a long name"; 
    }
    $ltype_name   = $stype2ltype_H{$stype_name};
    $lsource_name = $ssource2lsource_H{$ssource_name};
    if(! exists ($ltype_seqidx_H{$ltype_name})) { 
      die "ERROR type name $ltype_name not read in type sequence file, but read in map file, line: $line";
    }
    if(! exists ($lsource_seqidx_H{$lsource_name})) { 
      die "ERROR source name $lsource_name not read in source sequence file, but read in map file, line: $line";
    }
    if($mapstart < 0) { die "ERROR read non-positive start on line: $line"; }
    if($mapstop  < 0) { die "ERROR read non-positive stop on line: $line"; }
    my $type_len    = $ltype_len_H{$ltype_name};
    my $type_maplen = ($mapstart <= $mapstop) ? ($mapstop - $mapstart + 1) : ($mapstart - $mapstop + 1);
    my $source_len  = $lsource_len_H{$lsource_name};
    my $strand      = ($mapstart <= $mapstop) ? "+" : "-";

    $map_stype_idx_H{$stype_name} = $map_idx;
    $map_idx++;
    $map_stype_len_H{$stype_name} = $type_maplen;

    #if($mapstart > $source_len) { printf("WARNING: $line out of bounds $lsource_name length $source_len\n"); }
    #if($mapstop  > $source_len) { printf("WARNING: $line out of bounds $lsource_name length $source_len\n"); }
    #if($type_len != $type_maplen) { printf("WARNING: $line type sequence length coordinate length mismatch ($type_len != $type_maplen)\n"); }
    
    $stype_ssource_H{$stype_name} = $ssource_name;
    if(! exists ($ssource_stype_HA{$ssource_name})) { 
      @{$ssource_stype_HA{$ssource_name}} = ();
    }
    push(@{$ssource_stype_HA{$ssource_name}}, $stype_name);
    
    $stype_mapstart_H{$stype_name} = $mapstart;
    $stype_mapstop_H{$stype_name}  = $mapstop;
    $stype_strand_H{$stype_name}   = $strand;
#    if($strand ne "+") { 
#      printf("STRAND NEGATIVE $line");
#    }
  }
}
close(MAP);

# read the intron file
my %stype_intron_H = (); # key is short type sequence name, value is '1' if that sequence has an intron
read_intron_file($intron_file, \%stype_seqidx_H, \%stype_intron_H);

# extract the 16S sequences according to the map file from the source sequences 



# run cmsearch to predict 16S in the source fasta file
my $tblout_file   = $out_root . ".tblout";
my $cmsearch_file = $out_root . ".cmsearch";
my $cmsearch_cmd  = "cmsearch --rfam --cpu 0 --tblout $tblout_file $cm_file $source_file > $cmsearch_file";
#run_command($cmsearch_cmd, 0);

# remove overlaps 
my $deoverlap_tblout_file = $out_root . ".tblout.deoverlapped";
my $deoverlap_output_file = $out_root . ".tblout.deoverlapped.out";
my $deoverlap_cmd = "perl cmsearch-deoverlap.pl $tblout_file > $deoverlap_output_file";
#run_command($deoverlap_cmd, 1);

my $out_tbl_file               = $out_root . ".check.tbl";
my $sfetch_arc_file            = $out_root . ".check.RF01959.sfetch";
my $sfetch_bac_file            = $out_root . ".check.RF00177.sfetch";
my $sfetch_map_file            = $out_root . ".check.map.sfetch";
my $predicted_arc_fa_file      = $out_root . ".check.RF01959.predicted.fa";
my $predicted_bac_fa_file      = $out_root . ".check.RF00177.predicted.fa";
my $predicted_arc_stk_file     = $out_root . ".check.RF01959.predicted.stk";
my $predicted_bac_stk_file     = $out_root . ".check.RF00177.predicted.stk";
my $predicted_arc_cmalign_file = $out_root . ".check.RF01959.predicted.cmalign";
my $predicted_bac_cmalign_file = $out_root . ".check.RF00177.predicted.cmalign";
my $map_fa_file                = $out_root . ".check.map.fa";
my @out_tbl_AH = (); # array of hashes to fill with output information for the $out_tbl_file

# parse tblout and create output array and sfetch file
process_tblout_file($deoverlap_tblout_file, $sfetch_arc_file, $sfetch_bac_file, $sfetch_map_file, \%ssource_stype_HA, \%stype_ssource_H, \%lsource2ssource_H, \%stype2ltype_H, \%stype_mapstart_H, \%stype_mapstop_H, \%stype_strand_H, \%ssource_len_H, \%stype_len_H, \%map_stype_idx_H, \%map_stype_len_H, \%stype_intron_H, \@out_tbl_AH);

# fetch predicted sequences and use cmalign to align them
my $at_least_one_archaea  = 0;
my $at_least_one_bacteria = 0;
my $sfetch_cmd; 
my $cmalign_cmd;
#if(-e "$source_file.ssi") { run_command("rm $source_file.ssi", 1); }
#$sfetch_cmd = "esl-sfetch --index $source_file";
#run_command($sfetch_cmd, 1);

if(-s $sfetch_arc_file) { 
  $at_least_one_archaea = 1;
  $sfetch_cmd = "esl-sfetch -Cf $source_file $sfetch_arc_file > $predicted_arc_fa_file";
#  run_command($sfetch_cmd, 1);
  $cmalign_cmd = "cmfetch 2.cm RF01959 | cmalign --cpu 6 -o $predicted_arc_stk_file - $predicted_arc_fa_file > $predicted_arc_cmalign_file"; 
#  run_command($cmalign_cmd, 1);
}
if(-s $sfetch_bac_file) { 
  $at_least_one_bacteria = 1;
  $sfetch_cmd = "esl-sfetch -Cf $source_file $sfetch_bac_file > $predicted_bac_fa_file";
#  run_command($sfetch_cmd, 1);
  $cmalign_cmd = "cmfetch 2.cm RF00177 | cmalign --cpu 6 -o $predicted_bac_stk_file - $predicted_bac_fa_file > $predicted_bac_cmalign_file"; 
#  run_command($cmalign_cmd, 1);
}
# now fetch the sequences according to the map file coords:
$sfetch_cmd = "esl-sfetch -Cf $source_file $sfetch_map_file > $map_fa_file";
run_command($sfetch_cmd, 0);

# read the type fasta file into a sequence hash, and the newly created $map_fa_file into a sequence hash
my %type_seq_H = ();
my %map_seq_H  = ();
my $idx;
read_fasta($type_file, \%type_seq_H);
read_fasta($map_fa_file, \%map_seq_H);
#debug_print_hash(\%type_seq_H, "type_seq");
#debug_print_hash(\%map_seq_H, "map_seq");
# compare all sequences in both
foreach $stype_name (sort keys %map_stype_idx_H) { 
  $idx        = $map_stype_idx_H{$stype_name};
  $ltype_name = $stype2ltype_H{$stype_name};
  if(! exists $type_seq_H{$ltype_name})    { die "ERROR in fasta file, type sequence $ltype_name does not exist"; }
  if((! exists $map_seq_H{$ltype_name}) && 
     ($out_tbl_AH[($idx-1)]{"uf_MapLengthDiffersFromActual"} eq "") && 
     ($out_tbl_AH[($idx-1)]{"uf_MapCoordsOutOfBounds"} eq "") && 
     ($out_tbl_AH[($idx-1)]{"uf_NoInfernalPrediction"} eq "")) { 
    die "ERROR in fasta file, map sequence $ltype_name with identical length and no map errors does not exist"; 
  }
  if(exists $type_seq_H{$ltype_name} && exists $map_seq_H{$ltype_name}) { 
    if($type_seq_H{$ltype_name} ne $map_seq_H{$ltype_name}) { 
      $out_tbl_AH[($idx-1)]{"uf_MapSequenceDiffersFromActual"} = "VOID";
#      printf("$stype_name sequence differs");
    }  
    else { 
#      printf("$stype_name sequence identical");
    }
  }
}

# finally parse the cmalign file and get start and end coordinates of all predictions
if($at_least_one_archaea) { 
  parse_cmalign_file($predicted_arc_cmalign_file, \@out_tbl_AH, \%map_stype_idx_H);
}
if($at_least_one_bacteria) { 
  parse_cmalign_file($predicted_bac_cmalign_file, \@out_tbl_AH, \%map_stype_idx_H);
}
                     

# output tabular file
output_tabular_file($out_tbl_file, \@out_tbl_AH);
# 

#################################################################
# Subroutine : parse_seqstat_file()
# Incept:      EPN, Wed Dec 14 16:16:22 2016 [ribotyper-v1/ribo.pm]
#
# Purpose:     Parse an esl-seqstat -a output file.
#              
# Arguments: 
#   $seqstat_file:            file to parse
#   $nseq_R:                  REF to the number of sequences read, updated here
#   $seqidx_HR:               REF to hash of sequence indices to fill here
#   $len_HR:               REF to hash of sequence lengths to fill here
#
# Returns:     Total number of nucleotides read (summed length of all sequences). 
#              Fills %{$seqidx_HR} and %{$len_HR} and updates 
#              $$max_targetname_length_R, $$max_length_length_R, and $$nseq_R.
# 
# Dies:        If the sequence file has two sequences with identical names.
#              Error message will list all duplicates.
#              If no sequences were read.
#
################################################################# 
sub parse_seqstat_file { 
  my $nargs_expected = 4;
  my $sub_name = "parse_seqstat_file";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($seqstat_file, $nseq_R, $seqidx_HR, $len_HR) = @_;

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
      #printf("targetname: $targetname line: $line\n");
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
      $len_HR->{$targetname} = $length;
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

#################################################################
# Subroutine:  run_command()
# Incept:      EPN, Mon Dec 19 10:43:45 2016 [ribotyper-v1/ribo.pm]
#
# Purpose:     Runs a command using system() and exits in error 
#              if the command fails. If $be_verbose, outputs
#              the command to stdout. If $FH_HR->{"cmd"} is
#              defined, outputs command to that file handle.
#
# Arguments:
#   $cmd:         command to run, with a "system" command;
#   $be_verbose:  '1' to output command to stdout before we run it, '0' not to
#
# Returns:    amount of time the command took, in seconds
#
# Dies:       if $cmd fails
#################################################################
sub run_command {
  my $sub_name = "run_command()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmd, $be_verbose) = @_;
  
  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));

  system($cmd);

  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));

  if($? != 0) { 
    die "ERROR in $sub_name, the following command failed:\n$cmd\n";
  }

  return ($stop_time - $start_time);
}

#################################################################
# Subroutine : remove_dir_path()
# Incept:      EPN, Mon Nov  9 14:30:59 2009 [ssu-align] 
#
# Purpose:     Given a full path of a file remove the directory path.
#              For example: "foodir/foodir2/foo.stk" becomes "foo.stk".
#
# Arguments: 
#   $fullpath: name of original file
# 
# Returns:     The string $fullpath with dir path removed.
#
################################################################# 
sub remove_dir_path {
  my $sub_name = "remove_dir_path()";
  my $nargs_expected = 1;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my $fullpath = $_[0];

  $fullpath =~ s/^.+\///;

  return $fullpath;
}

#################################################################
# Subroutine : shorten_keys()
# Incept:      EPN, Mon Aug 14 11:13:57 2017
#
# Purpose:     Shorten keys in the idx and length hashes:
#              gi|756808612|ref|NG_041941.1| shortens to NG_041941
#
# Arguments: 
#   $idx_HR:  ref to hash, key is long sequence name, value is index in input file
#   $len_HR:  ref to hash: key is long sequence name, value is sequence length
#   $sidx_HR: ref to hash: key is short sequence name, value is index in input file
#   $slen_HR: ref to hash: key is short sequence name, value is sequence length
#   $s2l_HR:  ref to hash: key is short sequence name, value is long name
#   $l2s_HR:  ref to hash: key is long sequence name, value is short name
# 
# Returns:     void
# Dies:        if we try to make any duplicate keys
#
################################################################# 
sub shorten_keys { 
  my $sub_name = "shorten_keys()";
  my $nargs_expected = 6;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($idx_HR, $len_HR, $sidx_HR, $slen_HR, $s2l_HR, $l2s_HR) = @_;

  my $longkey;
  my $shortkey;
  foreach $longkey (keys %{$idx_HR}) { 
    if($longkey =~ m/^gi\|\d+\|\S+\|(\S+)\.\d+\|/) { 
      $shortkey = $1;
      $sidx_HR->{$shortkey} = $idx_HR->{$longkey};
      $l2s_HR->{$longkey} = $shortkey;
      if(exists $s2l_HR->{$shortkey}) { 
        die "ERROR in $sub_name, trying to store duplicate in s2l_HR->{$shortkey}, already exists as $s2l_HR->{$shortkey}";
      }
      $s2l_HR->{$shortkey} = $longkey;
      if(! exists $len_HR->{$longkey}) { 
        die "ERROR in $sub_name, no length value for $longkey"; 
      }
      $slen_HR->{$shortkey} = $len_HR->{$longkey};
    }
    else { 
      die "ERROR in $sub_name, key $longkey not in expected format";
    }
  }

  foreach my $longkey (keys %{$len_HR}) { 
    if(! exists $l2s_HR->{$longkey}) { 
      die "ERROR in $sub_name, long key $longkey only exists in idx"; 
    }
    if(! exists $slen_HR->{$l2s_HR->{$longkey}}) { 
      die "ERROR in $sub_name, short length not filled when reading idx"; 
    }
  }
  
  return;
}

#################################################################
# Subroutine : read_intron_file()
# Incept:      EPN, Mon Aug 14 15:51:31 2017
#
# Purpose:     Store info on sequences with known introns.
#
# Arguments: 
#   $intron_file:      file with list of accessions (short type strain sequence name)
#   $stype_idx_HR:     ref to hash, already filled with indexes
#   $stype_intron_HR:  ref to hash, key is short sequence name, value is 1 if 
#                      sequence has an intron in it.
# Returns:     void
#
################################################################# 
sub read_intron_file { 
  my $sub_name = "read_intron_file()";
  my $nargs_expected = 3;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($intron_file, $stype_idx_HR, $stype_intron_HR) = @_;

  open(IN, $intron_file) || die "ERROR, unable to open $intron_file";

  while(my $line = <IN>) { 
    chomp $line;
#    if(! exists $stype_idx_HR->{$line}) { 
#      die "ERROR in $sub_name, sequence $line does not exist in the short type sequence name index hash";
#    }
    $stype_intron_HR->{$line} = 1;
  }
  close(IN);

  return;
}

#################################################################
# Subroutine : process_tblout_file()
# Incept:      EPN, Tue Aug 15 09:23:40 2017
#
# Purpose:     Given a de-overlapped tblout file, create and print
#              the output
#
# Arguments: 
#   $tblout_file:        file with list of accessions (short type strain sequence name)
#   $sfetch_arc_file:    sfetch input file for archaeal sequences
#   $sfetch_bac_file:    sfetch input file for bacterial sequences
#   $sfetch_map_file:    sfetch input file for map sequences
#   $ssource_stype_HAR:  ref to hash of arrays, key is source short name, value is array of type short names
#                        that map to that source
#   $stype_ssource_HR:   ref to hash, key is short type name, value is short source name
#   $lsource2ssource_HR: ref to hash, key is long source name, value is short source name
#   $stype2ltype_HR:     ref to hash, key is short type name, value is long type name
#   $stype_mapstart_HR:  ref to hash, key is short type name, value is start position
#   $stype_mapstop_HR:   ref to hash, key is short type name, value is stop position
#   $stype_strand_HR:    ref to hash, key is short type name, value is strand
#   $stype_len_HR:       ref to hash, key is short type name, value is length
#   $ssource_len_HR:     ref to hash, key is short source name, value is length
#   $map_stype_idx_HR:   ref to hash, key is short type name, value is index in map file
#   $map_stype_len_HR:   ref to hash, key is short type name, value is length in map file
#   $stype_intron_HR:    ref to hash, key is short sequence name, value is 1 if 
#                        sequence has an intron in it.
#   $out_tbl_AHR:        ref to array of hashes to fill here with output lines for tabular output file
# Returns:     void
#
################################################################# 
sub process_tblout_file { 
  my $sub_name = "process_tblout_file()";
  my $nargs_expected = 17;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($tblout_file, $sfetch_arc_file, $sfetch_bac_file, $sfetch_map_file, $ssource_stype_HAR, $stype_ssource_HR, $lsource2ssource_HR, $stype2ltype_HR, $stype_mapstart_HR, $stype_mapstop_HR, $stype_strand_HR, $ssource_len_HR, $stype_len_HR, $map_stype_idx_HR, $map_stype_len_HR, $stype_intron_HR, $out_tbl_AHR) = @_;

  open(IN, $tblout_file) || die "ERROR, unable to open $tblout_file";

# for each line of tblout file, does this hit overlap with any of the sequences in
#@{$ssource_stype_HA{$source_name}}? (use $type_mapstart_H and $type_mapstop_H)
#If yes and we already output a line, die
#If yes and we haven't already output a line, output a line

  # open sfetch files
  open(ASFETCH, ">", $sfetch_arc_file) || die "ERROR unable to open $sfetch_arc_file for writing"; 
  open(BSFETCH, ">", $sfetch_bac_file) || die "ERROR unable to open $sfetch_bac_file for writing"; 
  open(MSFETCH, ">", $sfetch_map_file) || die "ERROR unable to open $sfetch_map_file for writing"; 

  my $uf_str = "";
  my $uf_noinfernalprediction_ct = 0;
  my $uf_maplengthdiffers_ct = 0;
  my $uf_mapstranddiffers_ct = 0;
  my $uf_predictiondisagrees_ct = 0;
  my $clean_ct = 0;
  my $stype_name;
  my $stype_class;
  my $nres_overlap;
  my $left_overhang;  
  my $right_overhang;
  my $domain;
  my $idx;
  my $has_intron;
  my $ssource_name;        
  my $do_output;           # flag telling us whether or not to output this sequence
  my $cur_mapstart;
  my $cur_mapstop;
  my $cur_strand;

  # values in a cmsearch tblout line
  my $seqfrom;
  my $seqto;
  my $target;
  my $model;
  my $acc;
  my $mdlfrom;
  my $mdlto;
  my $strand;
  my $score;
  my $evalue; 

  # initialize the output array of hashes, any sequence without a matching prediction will stay this way when output
  my %output_H = ();
  my $output_ct = 0;
  my @output_A = ();
  my $input_ct = scalar(keys %{$stype_len_HR});
  foreach $stype_name (sort keys %{$stype_len_HR}) { 
    if(! exists $map_stype_len_HR->{$stype_name}) { 
      die "ERROR stype_name $stype_name not in map_stype_len hash\n";
    }
    $idx = $map_stype_idx_HR->{$stype_name};
    $domain = "?";
    $has_intron = (exists $stype_intron_HR->{$stype_name}) ? "yes" : "no";
    $ssource_name = $stype_ssource_HR->{$stype_name};
    $stype_class = ($stype_mapstart_HR->{$stype_name} == 1) ? "full" : "sub"; 
    $cur_strand = $stype_strand_HR->{$stype_name};
    %{$out_tbl_AHR->[($idx-1)]} = ();
    $out_tbl_AHR->[($idx-1)]{"idx"}             = $idx;
    $out_tbl_AHR->[($idx-1)]{"domain"}          = "?";
    $out_tbl_AHR->[($idx-1)]{"stype_name"}      = $stype_name;
    $out_tbl_AHR->[($idx-1)]{"ssource_name"}    = $ssource_name;
    $out_tbl_AHR->[($idx-1)]{"stype_class"}     = $stype_class;
    $out_tbl_AHR->[($idx-1)]{"cur_strand"}      = $cur_strand;
    $out_tbl_AHR->[($idx-1)]{"pred_strand"}     = "?";
    $out_tbl_AHR->[($idx-1)]{"map_start"}       = $stype_mapstart_HR->{$stype_name};
    $out_tbl_AHR->[($idx-1)]{"map_stop"}        = $stype_mapstop_HR->{$stype_name};
    $out_tbl_AHR->[($idx-1)]{"pred_start"}      = "-";
    $out_tbl_AHR->[($idx-1)]{"pred_stop"}       = "-";
    $out_tbl_AHR->[($idx-1)]{"left_overhang"}   = "-";
    $out_tbl_AHR->[($idx-1)]{"right_overhang"}  = "-";
    $out_tbl_AHR->[($idx-1)]{"has_intron"}      = $has_intron;
    $out_tbl_AHR->[($idx-1)]{"pred_cmfrom"}     = "-";
    $out_tbl_AHR->[($idx-1)]{"pred_cmto"}       = "-";
    $out_tbl_AHR->[($idx-1)]{"uf_NoInfernalPrediction"}       = "VOID";
    $out_tbl_AHR->[($idx-1)]{"uf_PredictionDisagreesWithMap"}     = "";
    $out_tbl_AHR->[($idx-1)]{"uf_MapStrandDiffersFromPrediction"}  = "";
    $out_tbl_AHR->[($idx-1)]{"uf_MapSequenceDiffersFromActual"}   = "";
    if($stype_len_HR->{$stype_name} != $map_stype_len_HR->{$stype_name}) { 
      $out_tbl_AHR->[($idx-1)]{"uf_MapLengthDiffersFromActual"} = "(actual:" . $stype_len_HR->{$stype_name} . "!=map:" . $map_stype_len_HR->{$stype_name} . ");";
    }
    else { 
      $out_tbl_AHR->[($idx-1)]{"uf_MapLengthDiffersFromActual"} = "";
    }
    if(($stype_mapstart_HR->{$stype_name} > $ssource_len_HR->{$ssource_name}) ||
       ($stype_mapstop_HR->{$stype_name}  > $ssource_len_HR->{$ssource_name})) { 
      $out_tbl_AHR->[($idx-1)]{"uf_MapCoordsOutOfBounds"} = "(actual_source_len:" . $ssource_len_HR->{$ssource_name} . ",but_map:" . $stype_mapstart_HR->{$stype_name} . "-" . $stype_mapstop_HR->{$stype_name} . ");";
      #printf("$stype_name OUTOFBOUNDS\n");
    }      
    else { 
      $out_tbl_AHR->[($idx-1)]{"uf_MapCoordsOutOfBounds"}       = "";
    }
  }
  
  while(my $line = <IN>) { 
    chomp $line;
    if($line !~ m/^\#/) { 
      my @el_A = split(/\s+/, $line);

      ##target name             accession query name           accession mdl mdl from   mdl to seq from   seq to strand trunc pass   gc  bias  score   E-value inc description of target
      ##----------------------- --------- -------------------- --------- --- -------- -------- -------- -------- ------ ----- ---- ---- ----- ------ --------- --- ---------------------
      #lcl|dna_BP444_24.8k:251  -         SSU_rRNA_archaea     RF01959   hmm        3     1443        2     1436      +     -    6 0.53   6.0 1078.9         0 !   -
      if(scalar(@el_A) < 18) { die "ERROR found less than 18 columns in cmsearch tabular output at line: $line"; }
      ($target, $model, $acc, $mdlfrom, $mdlto, $seqfrom, $seqto, $strand, $score, $evalue) = 
          ($el_A[0], $el_A[2], $el_A[3], $el_A[5], $el_A[6], $el_A[7], $el_A[8], $el_A[9],  $el_A[14], $el_A[15]);

      if($score >= 50.) { 
        $ssource_name = $lsource2ssource_HR->{$target};
        if(! exists $ssource_stype_HA{$ssource_name}) { 
          die "ERROR couldn't find source $ssource_name from tblout out line: $line\n"; 
        }
        for(my $i = 0; $i < scalar(@{$ssource_stype_HAR->{$ssource_name}}); $i++) { 
          $stype_name = $ssource_stype_HA{$ssource_name}[$i];
          $stype_class = ($stype_mapstart_HR->{$stype_name} == 1) ? "full" : "sub"; 
          $nres_overlap = 0;
          # we don't enforce that strand must be identical
          if($strand eq "+") { 
            if($stype_strand_HR->{$stype_name} eq "+") { 
              $nres_overlap = get_overlap($seqfrom, $seqto, $stype_mapstart_HR->{$stype_name}, $stype_mapstop_HR->{$stype_name});
            }
            else { 
              $nres_overlap = get_overlap($seqfrom, $seqto, $stype_mapstop_HR->{$stype_name}, $stype_mapstart_HR->{$stype_name});
            }
          }
          else { # predicted strand is "-"
            if($stype_strand_HR->{$stype_name} eq "+") { 
              $nres_overlap = get_overlap($seqto, $seqfrom, $stype_mapstart_HR->{$stype_name}, $stype_mapstop_HR->{$stype_name});
            }
            else { 
              $nres_overlap = get_overlap($seqto, $seqfrom, $stype_mapstop_HR->{$stype_name}, $stype_mapstart_HR->{$stype_name});
            }
          }
          if($nres_overlap > 0) { 
            # sanity checks:
            $do_output = 1;
            if(exists $output_H{$stype_name}) { 
              printf"WARNING found two hits that overlap with $stype_name\n"; 
              $do_output = 0;
            }
            if(($acc ne "RF00177") && ($acc ne "RF01959")) { 
              die "ERROR unrecognized accession $acc"; 
            }
            if(! exists $map_stype_idx_HR->{$stype_name}) { 
              die "ERROR stype_name $stype_name not in map_stype_idx hash\n";
            }
            if(! exists $map_stype_len_HR->{$stype_name}) { 
              die "ERROR stype_name $stype_name not in map_stype_len hash\n";
            }
            if(! exists $stype_len_HR->{$stype_name}) { 
              die "ERROR stype_name $stype_name not in stype_len hash\n";
            }
            $cur_mapstart = ($strand eq $stype_strand_HR->{$stype_name}) ? $stype_mapstart_HR->{$stype_name} : $stype_mapstop_HR->{$stype_name};
            $cur_mapstop  = ($strand eq $stype_strand_HR->{$stype_name}) ? $stype_mapstop_HR->{$stype_name} : $stype_mapstart_HR->{$stype_name};
            $left_overhang  = $cur_mapstart - $seqfrom;
            $right_overhang = $cur_mapstop - $seqto;
            $domain = ($acc eq "RF00177") ? "bacteria" : "archaea";
            $idx = $map_stype_idx_HR->{$stype_name};
            $has_intron = (exists $stype_intron_HR->{$stype_name}) ? "yes" : "no";
            # check for errors:
            # 
            if($do_output) { 
              $out_tbl_AHR->[($idx-1)]{"idx"}             = $idx;
              $out_tbl_AHR->[($idx-1)]{"domain"}          = $domain;
              $out_tbl_AHR->[($idx-1)]{"stype_name"}      = $stype_name;
              $out_tbl_AHR->[($idx-1)]{"ssource_name"}    = $ssource_name;
              $out_tbl_AHR->[($idx-1)]{"stype_class"}     = $stype_class;
              $out_tbl_AHR->[($idx-1)]{"cur_strand"}      = $stype_strand_HR->{$stype_name};
              $out_tbl_AHR->[($idx-1)]{"pred_strand"}     = $strand;
              $out_tbl_AHR->[($idx-1)]{"map_start"}       = $cur_mapstart;
              $out_tbl_AHR->[($idx-1)]{"map_stop"}        = $cur_mapstop;
              $out_tbl_AHR->[($idx-1)]{"pred_start"}      = $seqfrom;
              $out_tbl_AHR->[($idx-1)]{"pred_stop"}       = $seqto;
              $out_tbl_AHR->[($idx-1)]{"left_overhang"}   = $left_overhang;
              $out_tbl_AHR->[($idx-1)]{"right_overhang"}  = $right_overhang;
              $out_tbl_AHR->[($idx-1)]{"has_intron"}      = $has_intron;
              $out_tbl_AHR->[($idx-1)]{"uf_NoInfernalPrediction"} = "";
              if($left_overhang != 0 || $right_overhang != 0) { 
                $out_tbl_AHR->[($idx-1)]{"uf_PredictionDisagreesWithMap"} = "(map:" . $cur_mapstart . "-" . $cur_mapstop . "!=prediction:" . $seqfrom . "-" . $seqto . ")"; 
              }
              if($stype_strand_HR->{$stype_name} ne $strand) { 
                $out_tbl_AHR->[($idx-1)]{"uf_MapStrandDiffersFromPrediction"} = "(predicted:" . $strand . "!=map:" . $stype_strand_HR->{$stype_name} . ");";
              }
              # already filled in values for uf_MapLengthDiffersFromActual and uf_MapCoordsOutOfBounds when initializing

              if($domain eq "archaea") { print ASFETCH ("$stype_name $seqfrom $seqto $target\n"); }
              else                     { print BSFETCH ("$stype_name $seqfrom $seqto $target\n"); }
              if(($out_tbl_AHR->[($idx-1)]{"uf_MapLengthDiffersFromActual"} eq "")  && 
                 ($out_tbl_AHR->[($idx-1)]{"uf_MapCoordsOutOfBounds"} eq "")) { 
                # only fetch sequences for which the length is identical and coordinates are sensical
                if(! exists $stype2ltype_HR->{$stype_name}) { 
                  die "ERROR in $sub_name, stype_name $stype_name not in stype2ltype_HR";
                }
                print MSFETCH ("$stype2ltype_HR->{$stype_name} $cur_mapstart $cur_mapstop $target\n");
              }
              #else { printf("SKIPPING $stype_name $cur_mapstart $cur_mapstop $target\n"); }
              $output_H{$stype_name} = 1;
            }
          }
        }
      }
    }
  }
  close(IN);
  close(ASFETCH);
  close(BSFETCH);
}

#################################################################
# Subroutine: get_overlap()
# Incept:     EPN, Mon Mar 14 13:47:57 2016 [dnaorg_scripts:dnaorg.pm:getOverlap()]
#
# Purpose:    Calculate number of nucleotides of overlap between
#             two regions.
#
# Args:
#  $start1: start position of hit 1 (must be <= $end1)
#  $end1:   end   position of hit 1 (must be >= $end1)
#  $start2: start position of hit 2 (must be <= $end2)
#  $end2:   end   position of hit 2 (must be >= $end2)
#
# Returns:  $noverlap:    Number of nucleotides of overlap between hit1 and hit2, 
#                         0 if none
#
# Dies:     if $end1 < $start1 or $end2 < $start2.
sub get_overlap {
  my $sub_name = "get_overlap";
  my $nargs_exp = 4;
  if(scalar(@_) != 4) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($start1, $end1, $start2, $end2) = @_; 

  # printf("in $sub_name $start1..$end1 $start2..$end2\n");

  if($start1 > $end1) { die "ERROR in $sub_name start1 > end1 ($start1 > $end1)"; }
  if($start2 > $end2) { die "ERROR in $sub_name start2 > end2 ($start2 > $end2)"; }

  # Given: $start1 <= $end1 and $start2 <= $end2.
  
  # Swap if nec so that $start1 <= $start2.
  if($start1 > $start2) { 
    my $tmp;
    $tmp   = $start1; $start1 = $start2; $start2 = $tmp;
    $tmp   =   $end1;   $end1 =   $end2;   $end2 = $tmp;
  }
  
  # 3 possible cases:
  # Case 1. $start1 <=   $end1 <  $start2 <=   $end2  Overlap is 0
  # Case 2. $start1 <= $start2 <=   $end1 <    $end2  
  # Case 3. $start1 <= $start2 <=   $end2 <=   $end1
  if($end1 < $start2) { return (0); }                    # case 1
  if($end1 <   $end2) { return ($end1 - $start2 + 1); }  # case 2
  if($end2 <=  $end1) { return ($end2 - $start2 + 1); }  # case 3
  die "ERROR in $sub_name, unforeseen case in $start1..$end1 and $start2..$end2";

  return; # NOT REACHED
}

#################################################################
# Subroutine : output_tabular_file()
# Incept:      EPN, Fri Aug 18 06:39:11 2017
#
# Purpose:     Output to the tabular output file, given an array 
#              of text to output.
#
# Arguments: 
#   $out_tbl_file:       file to output to
#   $out_tbl_AHR:         ref to array of lines to output
# Returns:     void
#
################################################################# 
sub output_tabular_file { 
  my $sub_name = "output_tabular_file()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($out_tbl_file, $out_tbl_AHR) = @_; 
  
  # output headers:
  open(OUT, ">", $out_tbl_file) || die "ERROR unable to open $out_tbl_file for writing"; 
  printf OUT ("%-5s  %-8s  %-10s  %-10s  %4s  %3s %3s  %21s   %21s    %21s  %11s  %6s  %18s\n", "#", "", "type", "", "full", "", "", "   current-coords", "  predicted-coords", "  coord-difference", "", "", "", "");
  printf OUT ("%-5s  %-8s  %-10s  %-10s  %4s  %3s %3s  %21s   %21s    %21s  %11s  %6s  %18s\n", "#", "pred", "strain", "source", "or", "cur", "prd", "---------------------", "---------------------", "---------------------", "predicted", "known", "");
  printf OUT ("%-5s  %-8s  %-10s  %-10s  %4s  %3s %3s  %10s %10s   %10s %10s    %10s %10s  %11s  %6s  %18s\n", "#idx", "domain", "accn", "accn", "sub", "str", "str", "start", "stop", "start", "stop", "start", "stop", "model-span", "intron", "unexpected_features");
  printf OUT ("%-5s  %-8s  %-10s  %-10s  %4s  %3s %3s  %10s %10s   %10s %10s    %10s %10s  %11s  %6s  %18s\n", "#----", "-------", "----------", "----------", "----", "---", "---", "----------", "----------", "----------", "----------", "----------", "----------", "-----------", "------", "------------------");

  my $i;
  my %ct_H = ();
  $ct_H{"input"}                         = scalar(@{$out_tbl_AHR});
  $ct_H{"predicted"}                     = 0;
  $ct_H{"uf_NoInfernalPrediction"}       = 0;
  $ct_H{"CLEAN"}                         = 0;
  $ct_H{"uf_PredictionDisagreesWithMap"} = 0;
  $ct_H{"uf_MapLengthDiffersFromActual"}    = 0;
  $ct_H{"uf_MapCoordsOutOfBounds"}      = 0;
  $ct_H{"uf_MapStrandDiffersFromPrediction"} = 0;
  $ct_H{"uf_MapSequenceDiffersFromActual"}  = 0;
  $ct_H{"CompleteModelSpan"}            = 0;

  for($idx = 1; $idx <= scalar(@{$out_tbl_AHR}); $idx++) { 
    my $uf_str     = create_ufeature_string(\%{$out_tbl_AHR->[($idx-1)]}, \%ct_H);
    my $model_span = create_modelspan_string(\%{$out_tbl_AHR->[($idx-1)]}, \%ct_H);
    printf OUT ("%-5s  %-8s  %-10s  %-10s  %4s  %3s %3s  %10s %10s   %10s %10s    %10s %10s  %11s  %6s  %s\n", 
                $out_tbl_AHR->[($idx-1)]{"idx"}, 
                $out_tbl_AHR->[($idx-1)]{"domain"}, 
                $out_tbl_AHR->[($idx-1)]{"stype_name"}, 
                $out_tbl_AHR->[($idx-1)]{"ssource_name"}, 
                $out_tbl_AHR->[($idx-1)]{"stype_class"}, 
                $out_tbl_AHR->[($idx-1)]{"cur_strand"}, 
                $out_tbl_AHR->[($idx-1)]{"pred_strand"}, 
                $out_tbl_AHR->[($idx-1)]{"map_start"}, 
                $out_tbl_AHR->[($idx-1)]{"map_stop"}, 
                $out_tbl_AHR->[($idx-1)]{"pred_start"}, 
                $out_tbl_AHR->[($idx-1)]{"pred_stop"}, 
                $out_tbl_AHR->[($idx-1)]{"left_overhang"},
                $out_tbl_AHR->[($idx-1)]{"right_overhang"},
                $model_span,
                $out_tbl_AHR->[($idx-1)]{"has_intron"},
                $uf_str);
  }
  close(OUT);
  $ct_H{"predicted"} = $ct_H{"input"} - $ct_H{"uf_NoInfernalPrediction"};

  print("\n# Saved tabular output to $out_tbl_file\n");

  print("\n# Summary statistics:\n");
  print("\n");
  printf("%-35s  %6s  %6s\n", "category", "count", "fraction");
  printf("%-35s  %6s  %6s\n", "--------------------", "------", "------");
  printf("%-35s  %6d  %6.4f\n", "Input",                     $ct_H{"input"},                        $ct_H{"input"}                        / $ct_H{"input"});
  printf("%-35s  %6d  %6.4f\n", "Predicted",                 $ct_H{"predicted"},                    $ct_H{"predicted"}                    / $ct_H{"input"});
  printf("%-35s  %6d  %6.4f\n", "NoInfernalPrediction",      $ct_H{"uf_NoInfernalPrediction"},      $ct_H{"uf_NoInfernalPrediction"}      / $ct_H{"input"});
  printf("%-35s  %6d  %6.4f\n", "CLEAN",                     $ct_H{"CLEAN"},                        $ct_H{"CLEAN"}                        / $ct_H{"input"});
  printf("%-35s  %6d  %6.4f\n", "PredictionDisagreesWithMap" ,   $ct_H{"uf_PredictionDisagreesWithMap"},    $ct_H{"uf_PredictionDisagreesWithMap"}    / $ct_H{"input"});
  printf("%-35s  %6d  %6.4f\n", "MapLengthDiffersFromActual",    $ct_H{"uf_MapLengthDiffersFromActual"},    $ct_H{"uf_MapLengthDiffersFromActual"}    / $ct_H{"input"});
  printf("%-35s  %6d  %6.4f\n", "MapSequenceDiffersFromActual",  $ct_H{"uf_MapSequenceDiffersFromActual"},  $ct_H{"uf_MapSequenceDiffersFromActual"}  / $ct_H{"input"});
  printf("%-35s  %6d  %6.4f\n", "MapCoordsOutOfBounds",      $ct_H{"uf_MapCoordsOutOfBounds"},      $ct_H{"uf_MapCoordsOutOfBounds"}      / $ct_H{"input"});
  printf("%-35s  %6d  %6.4f\n", "MapStrandDiffersFromPrediction", $ct_H{"uf_MapStrandDiffersFromPrediction"}, $ct_H{"uf_MapStrandDiffersFromPrediction"} / $ct_H{"input"});
  printf("\n");
  printf("%-35s  %6d  %6.4f\n", "PredictionCompleteModelSpan", $ct_H{"CompleteModelSpan"}, $ct_H{"CompleteModelSpan"} / $ct_H{"input"});
  return;
}

#################################################################
# subroutine : read_fasta
# sub class  : crw and sequence
# 
# EPN 03.08.05 [Janelia: M_seq.pm]
#
# purpose : Open, read, and store the information in a given
#           .fa (fasta format) file.
#
# args : (1) $in_file
#            name of .fa file in current directory
#        (2) $seq_hash_ref
#            reference to the hash that will contain the sequence
#            information.  Fasta description line used as key for
#            each sequence, sequence is value.
################################################################# 
sub read_fasta
{
  my $sub_name = "read_fasta()";
  my $nargs_expected = 2;
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 
  my ($in_file, $seq_hash_ref) = @_;
  if(! -e $in_file) { die "ERROR in $sub_name, $in_file does not exist"; }
  if(! -s $in_file) { die "ERROR in $sub_name, $in_file exists but is empty"; }
  open(IN, $in_file) || die "ERROR unable to open $in_file";
    
  my $seqname = undef;
  my $line;

  while($line = <IN>) { 
    chomp $line;
    if($line ne "") { 
      if($line =~ /^\>(\S+)/) { 
        $seqname = $1;
        $seq_hash_ref->{$seqname} = "";
      }
      else { 
        if(defined $seqname) { 
          $seq_hash_ref->{$seqname} .= $line; 
        }
      }
    }
  }
  return;
}

#################################################################
# Subroutine: create_ufeature_string()
# Incept:     EPN, Fri Aug 18 10:01:25 2017
#
# Purpose:    Given an %output_H, create the ufeature string.
#
# Args:
#  %output_HR: hash of output values for current sequence
#  $ct_HR:     ref to count hash of unexpected features
#
# Returns:  $uf_str: string describing unexpected features
#
sub create_ufeature_string { 
  my $sub_name = "create_ufeature_string";
  my $nargs_exp = 2;
  if(scalar(@_) != 2) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($output_HR, $ct_HR) = @_; 

  my $uf_str = "";
  foreach my $uf ("uf_NoInfernalPrediction", 
                  "uf_PredictionDisagreesWithMap",
                  "uf_MapLengthDiffersFromActual",
                  "uf_MapSequenceDiffersFromActual",
                  "uf_MapStrandDiffersFromPrediction",
                  "uf_MapCoordsOutOfBounds") { 
    if(! exists $output_HR->{$uf}) { 
      die "ERROR in $sub_name, output_HR $uf does not exist"; 
    }
    if(! exists $ct_HR->{$uf}) { 
      die "ERROR in $sub_name, ct_HR $uf does not exist"; 
    }
    if($output_HR->{$uf} ne "") { 
      my $uf_toprint = $uf;
      $uf_toprint =~ s/^uf\_//;
      $uf_str .= $uf_toprint;
      if($output_HR->{$uf} ne "VOID") { 
        $uf_str .= $output_HR->{$uf};
      }
      $uf_str .= ";";
      $ct_HR->{$uf}++;
    }
  }

  if($uf_str eq "") { 
    $ct_HR->{"CLEAN"}++;
    $uf_str = "-";
  }
  return $uf_str;
}

#################################################################
# Subroutine: create_modelspan_string()
# Incept:     EPN, Fri Aug 18 11:42:29 2017
#
# Purpose:    Given an %output_H, create the modelspan string.
#
# Args:
#  %output_HR: hash of output values for current sequence
#  $ct_HR:     ref to count hash of unexpected features
#
# Returns:  $modelspan_str: string describing modelspan
#
sub create_modelspan_string { 
  my $sub_name = "create_modelspan_string";
  my $nargs_exp = 2;
  if(scalar(@_) != 2) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($output_HR, $ct_HR) = @_; 

  my $modelspan_str = "";
  if($output_HR->{"domain"} eq "archaea") { 
    $modelspan_str =  $output_HR->{"pred_cmfrom"} . "-" . $output_HR->{"pred_cmto"};
    if($output_HR->{"pred_cmfrom"} eq "-" && $output_HR->{"pred_cmto"} eq "-") { 
      $modelspan_str = "-";
    }
    elsif($output_HR->{"pred_cmfrom"} == 1 && $output_HR->{"pred_cmto"} == 1477) { 
      $modelspan_str .= "[]";
      $ct_HR->{"CompleteModelSpan"}++;
    }
    elsif($output_HR->{"pred_cmfrom"} == 1) { 
      $modelspan_str .= "[.";
    }
    elsif($output_HR->{"pred_cmto"} == 1477) { 
      $modelspan_str .= ".]";
    }
    else { 
      $modelspan_str .= "..";
    }
  }
  elsif($output_HR->{"domain"} eq "bacteria") { 
    $modelspan_str =  $output_HR->{"pred_cmfrom"} . "-" . $output_HR->{"pred_cmto"};

    if($output_HR->{"pred_cmfrom"} eq "-" && $output_HR->{"pred_cmto"} eq "-") { 
      $modelspan_str = "-";
    }
    elsif($output_HR->{"pred_cmfrom"} == 1 && $output_HR->{"pred_cmto"} == 1533) { 
      $modelspan_str .= "[]";
      $ct_HR->{"CompleteModelSpan"}++;
    }
    elsif($output_HR->{"pred_cmfrom"} == 1) { 
      $modelspan_str .= "[.";
    }
    elsif($output_HR->{"pred_cmto"} == 1533) { 
      $modelspan_str .= ".]";
    }
    else { 
      $modelspan_str .= "..";
    }
  }
  else { 
    $modelspan_str = "-";
  }
  
  return $modelspan_str;
}

#################################################################
# subroutine : debug_print_hash
# sub class  : general
# 
# EPN 03.08.05
# 
# purpose : Print to standard output the keys and values of a 
#           given hash
#
# args : (1) $hash_ref 
#            reference to hash to print
#        (2) $hash_name
#            name of hash to print
################################################################# 

sub debug_print_hash
{
  my ($hash_ref, $hash_name) = @_;
    
  print("IN DEBUG PRINT HASH\n");
  print("printing hash : $hash_name\n");
  my $i = 1;
  foreach my $header (sort keys (%{$hash_ref}))
  {
    print("$i KEY    : $header\n");
    print("$i VALUE : $hash_ref->{$header}\n");
    $i++;
  }
  print("finished printing hash : $hash_name\n");
  print("LEAVING DEBUG PRINT HASH\n");
}

#################################################################
# Subroutine : parse_cmalign_file()
# Incept:      EPN, Fri Aug 18 11:32:24 2017
#
# Purpose:     Parse a cmalign file, storing only the start and end 
#              positions in the model in @{$out_tbl_AHR}
# 
#              
# Arguments: 
#   $cmalign_file: file to parse
#   $out_tbl_AHR:  ref to array of hashes, with output info
#   $stype_idx_HR: ref to hash, key is short type name, value is index in map file
#
# Returns:     void; 
#
################################################################# 
sub parse_cmalign_file { 
  my $nargs_expected = 3;
  my $sub_name = "parse_cmalign_file";
  if(scalar(@_) != $nargs_expected) { printf STDERR ("ERROR, $sub_name entered with %d != %d input arguments.\n", scalar(@_), $nargs_expected); exit(1); } 

  my ($cmalign_file, $out_tbl_AHR, $map_stype_idx_HR) = @_;

  open(IN, $cmalign_file) || die "ERROR unable to open cmalign file $cmalign_file for reading";

  while(my $line = <IN>) { 
##                                                                           running time (s)                 
##                                                                    -------------------------------          
## idx  seq name   length  cm from    cm to  trunc    bit sc  avg pp  band calc  alignment      total  mem (Mb)
## ---  ---------  ------  -------  -------  -----  --------  ------  ---------  ---------  ---------  --------
#    1  NR_043409    1493        1     1477     no   1501.40   0.987       0.37       0.19       0.56     50.52
#    2  NR_043410    1497        1     1477     no   1541.35   0.989       0.36       0.19       0.55     50.55
#    3  NR_029127    1496        1     1477     no   1568.97   0.987       0.38       0.17       0.55     50.16
    chomp $line; 
    if($line !~ /^\#/) { 
      $line =~ s/^\s+//; # remove leading whitespace
      $line =~ s/\s+$//; # remove trailing whitespace
      my @el_A = split(/\s+/, $line);
      if(scalar(@el_A) != 12) { die "ERROR in $sub_name, unexpected number of tokens on cmalign output file line: $line";  }
      my ($stype_name, $cmfrom, $cmto) = ($el_A[1], $el_A[3], $el_A[4]);
      if(! exists $map_stype_idx_HR->{$stype_name}) { 
        die "ERROR in $sub_name, $stype_name does not exist in stype_idx_HR";
      }
      my $idx = $map_stype_idx_HR->{$stype_name};
      $out_tbl_AHR->[($idx-1)]{"pred_cmfrom"} = $cmfrom;
      $out_tbl_AHR->[($idx-1)]{"pred_cmto"}   = $cmto;
    }
  }
  close(IN);
  return;
}
