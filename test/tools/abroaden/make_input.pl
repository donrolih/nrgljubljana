#!/usr/bin/env perl
# Write test spectra for abroaden and broaden: a fixed set of delta peaks of total weight 1 at frequencies spread over
# several decades, both signs.
#
#   make_input.pl abin FILE [SIGMA]       records (omega, sigma, w); SIGMA (default ln 2) for every peak,
#                                         or "spread" for a sigma that varies from peak to peak
#   make_input.pl bin FILE                records (omega, w), the same peaks, for broaden
#   make_input.pl abin FILE SIGMA SCALE   weights multiplied by SCALE
use strict;
use warnings;

my ($kind, $file, $sigma, $scale) = @ARGV;
die "usage: $0 abin|bin FILE [SIGMA|spread] [SCALE]\n" unless defined $file && ($kind eq "abin" || $kind eq "bin");
$sigma //= log(2);
$scale //= 1.0;
my @peaks;
my $n = 0;
for my $e (1e-5, 3e-5, 1e-4, 4e-4, 1e-3, 3e-3, 1e-2, 2.5e-2, 0.07, 0.2, 0.45, 0.9) {
  for my $s (1, -1) {
    my $w = 1.0 + 0.1 * ($n % 7);
    push @peaks, [$s * $e, $w];
    $n++;
  }
}
my $total = 0;
$total += $_->[1] for @peaks;
open(my $f, ">:raw", $file) or die "cannot write $file: $!\n";
$n = 0;
for my $p (@peaks) {
  my ($e, $w) = ($p->[0], $scale * $p->[1] / $total);
  if ($kind eq "bin") {
    print $f pack("d2", $e, $w);
  } else {
    my $s = $sigma eq "spread" ? 0.1 + 0.05 * ($n % 13) : $sigma;
    print $f pack("d3", $e, $s, $w);
  }
  $n++;
}
close $f;
