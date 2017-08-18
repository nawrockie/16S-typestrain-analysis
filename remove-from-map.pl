use strict;

my $usage = "perl remove-from-map.pl <file with list of accession to remove> <map file>";
if(scalar(@ARGV) != 2) { die $usage; }

my ($remove_file, $map_file) = (@ARGV);

open(REMOVE, $remove_file) || die "ERROR, unable to open $remove_file";
my $line;
my $name;
my %remove_H = ();
while($line = <REMOVE>) { 
  chomp $line;
  if($line !~ m/^\#/) { 
    $line =~ s/\s+.*$//; # remove everything after first space, including space
    if($line =~ /gi\|\d+\|\S+\|(\S+)\.\d+|\S*$/) { 
      $name = $1;
      $remove_H{$name} = 1;
#      printf("removing $name\n");
    }
    else { 
      die "ERROR, couldn't parse $line";
    }
  }
}
close(REMOVE);

open(MAP, $map_file) || die "ERROR unable to open $map_file";
while($line = <MAP>) { 
  chomp $line;
  if($line !~ m/^\#/) { 
    if($line =~ /^(N\S+)\|\S+\.\d+\:\d+\-\d+\|/) { 
      $name = $1;
    }
    elsif($line =~ /^(N\S+)\|\S+\:\d+\-\d+\|/) { 
      $name = $1;
    }
    else { 
      die "ERROR unable to parse map line: $line";
    }
    if(! exists $remove_H{$name}) { 
      print $line . "\n";
    }
  }
}
