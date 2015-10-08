#!/usr/bin/perl -w

my %effects;

my $univariate_analysis = shift;
my $trivariate_analysis = shift;
my $th = shift || 0.05;

open(UNI, $univariate_analysis);
open(TRI, $trivariate_analysis);

my %pvals;

while(<UNI>) {
  chomp;
  my($enh, $pos, $effect, $p) = split(/\t/, $_);
  $pvals{$enh}{$pos} = $p;
}


while(<TRI>) {
  chomp;
  my($enh,$pos,$mut,$effect,$p) = split(/\t/, $_);
  $counts{$enh}{$pos}{$mut} = $effect;
}

print "enhancer\tposition\tmutation\teffect\ttstv\tuni_p\n";

while(my($enh, $pos_hash_ref) = each %counts) {
  my %pos_hash = %{$pos_hash_ref};
  while(my($pos, $mut_hash_ref) = each %pos_hash) {
    my %pos_hash = %{$pos_hash_ref};

    if(!exists($counts{$enh}{$pos}{A})) { $counts{$enh}{$pos}{A} = 0; }
    if(!exists($counts{$enh}{$pos}{C})) { $counts{$enh}{$pos}{C} = 0; }
    if(!exists($counts{$enh}{$pos}{G})) { $counts{$enh}{$pos}{G} = 0; }
    if(!exists($counts{$enh}{$pos}{T})) { $counts{$enh}{$pos}{T} = 0; }

    $effects{$enh}{$pos}{AC} = abs($counts{$enh}{$pos}{A} - $counts{$enh}{$pos}{C});
    $effects{$enh}{$pos}{AG} = abs($counts{$enh}{$pos}{A} - $counts{$enh}{$pos}{G});
    $effects{$enh}{$pos}{AT} = abs($counts{$enh}{$pos}{A} - $counts{$enh}{$pos}{T});
    $effects{$enh}{$pos}{CG} = abs($counts{$enh}{$pos}{C} - $counts{$enh}{$pos}{G});
    $effects{$enh}{$pos}{CT} = abs($counts{$enh}{$pos}{C} - $counts{$enh}{$pos}{T});
    $effects{$enh}{$pos}{GT} = abs($counts{$enh}{$pos}{G} - $counts{$enh}{$pos}{T});

    print "$enh\t$pos\tAC\t$effects{$enh}{$pos}{AC}\tTv\t$pvals{$enh}{$pos}\n";
    print "$enh\t$pos\tAG\t$effects{$enh}{$pos}{AG}\tTs\t$pvals{$enh}{$pos}\n";
    print "$enh\t$pos\tAT\t$effects{$enh}{$pos}{AT}\tTv\t$pvals{$enh}{$pos}\n";
    print "$enh\t$pos\tCG\t$effects{$enh}{$pos}{CG}\tTv\t$pvals{$enh}{$pos}\n";
    print "$enh\t$pos\tCT\t$effects{$enh}{$pos}{CT}\tTs\t$pvals{$enh}{$pos}\n";
    print "$enh\t$pos\tGT\t$effects{$enh}{$pos}{GT}\tTv\t$pvals{$enh}{$pos}\n";
  }
}
