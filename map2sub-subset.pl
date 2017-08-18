while($line = <>) { 
  if($line !~ m/^\#/) { 
    chomp $line;
    if($line =~ /^N\S+\|\S+\:(\d+)\-\d+\|/) { 
      if($1 ne "1") { 
        print $line . "\n";
      }
    }
    else { 
      die "ERROR $line is unexpected format";
#      printf("$line\n");
    }
  }
}
