$usage = "perl prune-tblout.pl <list of accessions to keep> <tblout file>";
if(scalar(@ARGV) != 2) { die $usage; }

($list_file, $tblout_file) = (@ARGV);
%keep_H = ();
open(IN, $list_file) || die "ERROR unable to open list file $list_file";
while($line = <IN>) { 
  chomp $line;
  $keep_H{$line} = 1;
}
close(IN);

open(TBLOUT, $tblout_file) || die "ERROR unable to open tblout file $tblout_file"; 
while($line = <TBLOUT>) { 
  $orig_line = $line;
  chomp $line;
  if($line !~ m/^\#/) { 
    $line =~ s/\s+.*$//;
    if(exists $keep_H{$line}) { 
      print $orig_line;
    }
  }
  else { 
    print $orig_line;
  }
}
close(TBLOUT);
