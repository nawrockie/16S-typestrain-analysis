while($line = <>) { 
  if($line !~ m/^\#/) { 
    chomp $line;
    if($line =~ /^N\S+\|(\S+)\:\d+\-\d+\|/) { 
      print $1 . "\n";
    }
    else { 
      die "ERROR $line is unexpected format";
#      printf("$line\n");
    }
  }
}
