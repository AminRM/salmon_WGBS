#!/usr//bin/perl -w

my $method = "match";
my $usage = "USAGE: $0 <gtf file> <snp file>\n";
my $reffile = shift or die $usage;
my $snpfile = shift or die $usage;
#my $goup = 2000;
#my $godown = 500;

my %genehash = ();


$reffile = "gunzip -c $reffile |" if ($reffile =~ /\.gz$/);
open(GE,$reffile) or die "Cannot open $reffile\n";
while(<GE>){
  next if(/^#/);
  chomp;
  my $loc;
  #my($acc,$chr,$gstart,$gend,$frm,$gstrand,$scr,$tend,$exons,$exstart,$exend,$src,$feat,$frame,$score,$rest)  = split(/\t/,$_);
  #my($chr,$src,$feat,$gstart,$gend,$scr,$gstrand,$frm,$rest) = split(/\t/,$_);
  my($taxid,$org,$gid,$cid,$status,$gene,$alias,$desc,$des,$maploc,$chr,$acc,$gstart,$gend,$gstrand,$exons,$omim) = split(/\t/,$_);
  if ($gstrand eq "plus"){$gstrand = "+";}
  elsif ($gstrand eq "minus"){$gstrand = "-";}
 
  my $rec = {
             GSTART => $gstart,
             GEND => $gend,
             STRAND => $gstrand,             
             GENE => $gene,
             GID => $gid,
             DESC => $desc,
             ACC => $acc,
             };

  my $coords = $gstart."_".$gend;
  $genehash{$chr}{$coords} = $rec; 
# print "GENE:$gene\t$gid\t$chr\t$gstart\t$gend\t$gstrand\n";
}
close(GE);
#exit;
$snpfile = "gunzip -c $snpfile |" if ($snpfile =~ /\.gz$/);
open(IN,$snpfile) or die "Cannot open $snpfile\n";
print "SNP\tCHR\tSTART\tEND\tDIST\tGENE\tGENE_START\tGENE_END\tSTRAND\tGENEID\tLOC\n";
while(<IN>){
   next if(/^#/);
   next if(/SNP/);
  chomp;
  my $dist = 0;
  my $closestDist = 9999999999;
  my $nextDist = $closestDist;
  my $closest = '.';
  my $nextClosest = '.';
  my $closestLoc = '.';
  my ($chr,$start,$stop,$rest) = split(/\t/,$_);
#print "SNP:$chr\t$start\t$stop\n";
  if ($chr eq "0"){
     print "$chr\t$start\t$stop\t.\t.\t.\t.\t$closest\t$closestLoc\n";
     next;
  }

#warn "$acc\t$chr\t$coord\n";
  my($genestart,$geneend,$gstrand,$gfeat,$gtrans);
#  next if ($chr =~ /^chrUn/);
  foreach my $key (keys %{$genehash{$chr}}){
      #print "$chr\t$coord\t$key\n";
      my ($tstart,$tend) = split(/_/,$key);
      $genestart = $genehash{$chr}{$key}->{GSTART};
      $geneend = $genehash{$chr}{$key}->{GEND};
      $gstrand = $genehash{$chr}{$key}->{STRAND};
#      $geneid = $genehash{$chr}{$key}->{GID};
#      $gfeat = $genehash{$chr}{$key}->{FEAT};
#      $gtrans = $genehash{$chr}{$key}->{TRANSCRIPT};
#      next if ($gfeat eq "gene");
      #$readNum = $pos_read_num - $neg_read_num;
      if($gstrand eq '+'){
         if($start < $tstart && $stop < $tstart){
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
      if ($dist < $closestDist){
             $nextDist= $closestDist;
             $closestDist = $dist;
             $closestLoc = $loc;
             #$closest = "$genehash{$chr}{$key}->{ACC}\t$chr\t$tstart\t$tend\t$gstrand\t$genehash{$chr}{$key}->{GENE}";
             $nextClosest = $closest;
             #$closest = "$genehash{$chr}{$key}->{GENE}\t$genehash{$chr}{$key}->{ACC}\t$tstart\t$tend\t$gstrand\t$genehash{$chr}{$key}->{FEAT}";
             $closest = "$genehash{$chr}{$key}->{GENE}\t$genehash{$chr}{$key}->{GSTART}\t$genehash{$chr}{$key}->{GEND}\t$gstrand\t$genehash{$chr}{$key}->{ACC}\t$genehash{$chr}{$key}->{DESC}";
      }
      elsif ($dist < $nextDist){
            $nextClosest = "$genehash{$chr}{$key}->{GENE}\t$genehash{$chr}{$key}->{GID}\t$gstrand";
      }
   }
   if($closestLoc eq ''){
      $closestDist = 100000000;
      $closest = "\t\t\t";
   }
#   print "$closest\t$genestart\t$geneend\t$coord\t$closestDist\t$closestLoc\t$rank\n";
#   print "$rank\t$name\t$chrn\t$coord\t$all1\t$all2\t$closestDist\t$closest\t$closestLoc\n";
#   print "$acc\t$chr\t$coord\t$closestDist\t$closest\t$closestLoc\t$nextClosest\t$nextDist\n";
   print "$chr\t$start\t$stop\t$closestDist\t$closest\t$closestLoc\n";
}
close(IN);

