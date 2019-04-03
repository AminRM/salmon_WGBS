#!/usr//bin/perl -w

my $method = "match";
my $usage = "USAGE: $0 <gtf file> <snp file>\n";
my $reffile = shift or die $usage;
my $snpfile = shift or die $usage;
#my $goup = 2000;
#my $godown = 500;

my %genehash = ();
my %exonhash = ();

$reffile = "gunzip -c $reffile |" if ($reffile =~ /\.gz$/);
open(GE,$reffile) or die "Cannot open $reffile\n";
while(<GE>){
  my($acc,$chr,$gstart,$gend,$frm,$gstrand,$scr,$tend,$exon,$exstart,$exend,$src,$feat,$frame,$score,$rest);
  next if(/^#/);
  chomp;
  my $loc;
  #my($taxid,$org,$gid,$cid,$status,$gene,$alias,$desc,$des,$maploc,$chr,$acc,$gstart,$gend,$gstrand,$exons,$omim) = split(/\t/,$_);
  ($chr,$src,$feat,$gstart,$gend,$scr,$gstrand,$frm,$rest) = split(/\t/,$_);
  if ($gstrand eq "plus"){$gstrand = "+";}
  elsif ($gstrand eq "minus"){$gstrand = "-";}
  next if ($feat =~ /transcript/);
  next unless ($feat =~ /exon/ || $feat =~ /mRNA/);

  if ($rest =~ /Genbank:(.*?)\;/) {
     $gid = $1;
    }
  #if ($rest =~ /transcript_id\=(.*?)/) {
  #   $gid = $1;
  #  }

  
  if ($rest =~ /ID=(.*?)\;/) {
     $exon = $1;
    }
 if ($feat =~ /mRNA/){
   my $grec = {
             GSTART => $gstart,
             GEND => $gend,
             STRAND => $gstrand,             
             GID => $gid,
             };

   my $coords = $gstart."_".$gend;
   $genehash{$chr}{$coords} = $grec; 
  # warn "GENE GID$gid\t$chr\t$gstart\t$gend\t$gstrand\n";
  }
  if ($feat =~ /exon/){
      my $xrec = {
              XSTART => $gstart,
              XEND => $gend,
              XSTRAND => $gstrand,
              XGID => $gid,
              };
    my $xcoords = $gstart."_".$gend;
    $exonhash{$gid}{$xcoords} = $xrec;
  }
}
close(GE);
#exit;
$snpfile = "gunzip -c $snpfile |" if ($snpfile =~ /\.gz$/);
open(IN,$snpfile) or die "Cannot open $snpfile\n";
#print "SNP\tCHR\tCOORD\tDIST\tGENE\tGENE_START\tGENE_END\tSTRAND\tGENEID\tLOC\n";
while(<IN>){
   next if(/^#/);
   next if(/ADD/);
  chomp;
  my $closestDist = 9999999999;
  my $nextDist = $closestDist;
  my $closest = '';
  my $closestgid = '';
  my $nextClosest = '';
  my $closestLoc = '';
#  my ($snp,$gchr,$coord,$gdist,$gene,$chr,$start,$end,$strand,$gid,$loc,$name) = split(/\t/,$_);
  my ($chr,$start,$stop,$rest) = split(/\t/,$_);


  my($genestart,$geneend,$gstrand,$gfeat,$gtrans);
  foreach my $key (keys %{$genehash{$chr}}){
      my $dist = 0;
      my ($tstart,$tend) = split(/_/,$key);
      $genestart = $genehash{$chr}{$key}->{GSTART};
      $geneend = $genehash{$chr}{$key}->{GEND};
      $gstrand = $genehash{$chr}{$key}->{STRAND};
      $gid = $genehash{$chr}{$key}->{GID};
    if($gstrand eq '+'){
         if($start< $tstart && $stop < $tstart){
             $dist = $tstart - $stop;
             $loc = "downstream";
         }
         elsif($start > $tend && $stop > $tend){
             $dist = $start - $tend;
             $loc = "upstream";
         }
         else{
              $dist = 0;
              $loc = "overlap";
        }
      }
      else{
         if($start < $tstart && $stop < $tstart){
             $dist = $tstart - $start;
             $loc = "upstream";
         }
         elsif($start > $tend && $stop > $tend){
             $dist = $stop - $tend;
             $loc = "downstream";
         }
          else{
              $dist = 0;
              $loc = "overlap";
        }

     }
     #print "DIST $gid  $dist\n";
      if ($dist < $closestDist){
             $nextDist= $closestDist;
             $closestDist = $dist;
             $closestLoc = $loc;
             $nextClosest = $closest;
             $closest = "$genehash{$chr}{$key}->{GID}\t$genehash{$chr}{$key}->{GSTART}\t$genehash{$chr}{$key}->{GEND}\t$gstrand";
             $closestgid = $genehash{$chr}{$key}->{GID};
      }
      elsif ($dist < $nextDist){
            $nextClosest = "$genehash{$chr}{$key}->{GID}\t$gstrand";
      }
   }
   if($closestLoc eq ''){
      $closestDist = 10000000000;
      $closest = "\t\t\t";
   }

 #  print "$chr\t$start\t$stop\t$closestDist\t$closest\t$closestLoc\n";
  if ($closestLoc eq "overlap"){
      my $xclosestDist = 99999999999;
      foreach my $key (keys %{$exonhash{$closestgid}}){
      my $xdist = 0;
      #print "$chr\t$coord\t$key\n";
      my ($tstart,$tend) = split(/_/,$key);
      $tstart = $exonhash{$closestgid}{$key}->{XSTART};
      $tend = $exonhash{$closestgid}{$key}->{XEND};
      $xstrand = $exonhash{$closestgid}{$key}->{XSTRAND};
    if($xstrand eq '+'){
         if($start< $tstart && $stop < $tstart){
             $dist = $tstart - $stop;
             $loc = "intron";
         }
         elsif($start > $tend && $stop > $tend){
             $dist = $start - $tend;
             $loc = "intron";
         }
         else{
              $dist = 0;
              $loc = "exon";
        }
      }
      else{
         if($start < $tstart && $stop < $tstart){
             $dist = $tstart - $start;
             $loc = "intron";
         }
         elsif($start > $tend && $stop > $tend){
             $dist = $stop - $tend;
             $loc = "intron";
         }
          else{
              $dist = 0;
              $loc = "exon";
        }

     }

      if ($xdist < $xclosestDist){
             $xclosestDist = $xdist;
             $closestLoc = $loc;
             #$closest = "$exonhash{$closestgid}{$key}->{XSTART}\t$exonhash{$closestgid}{$key}->{XEND}";
             $closest = "$exonhash{$closestgid}{$key}->{XGID}\t$exonhash{$closestgid}{$key}->{XSTART}\t$exonhash{$closestgid}{$key}->{XEND}\t$xstrand";
      }
      }
   }

  
   print "$chr\t$start\t$stop\t$closestDist\t$closest\t$closestLoc\n";
}
close(IN);

