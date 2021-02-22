system("wget ftp://hgdownload.cse.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz")

my %chr;
open(chrname,"name_chr.txt");
while ( <chrname> ){
  chomp $_;
  $chr{$_} = 1;
}
close (chrname);

my $file = "rheMac10.fa.gz";
open(in, "zcat $file |");
while ( <in> ) {
  chomp $_;
  if ($_ =~ />/) {
    $name = $_;
    $ckpt = 0;
    $name =~ s/>//g;
    if ($chr{$name} ne "") {
      print ($name,"\n");
      $ckpt = 1;
      open($out,"|gzip -c > $name.fa.gz");
      print $out (">$name\n");
    } else {
      close ($out);
    }
  } else {
    if ($ckpt == 1) {
      print $out ("$_\n");   
    }
  }   
}
close(in);
