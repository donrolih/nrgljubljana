#!/usr/bin/env perl
# Compare two spec.dat files (two columns, the same mesh): |y1 - y2| <= TOL * |y1| + FLOOR at every point.
#
#   compare_spec.pl --tol=TOL --floor=FLOOR file1 file2
use strict;
use warnings;
use Getopt::Long;

my ($tol, $floor) = (1e-10, 1e-14);
GetOptions("tol=f" => \$tol, "floor=f" => \$floor) or die "bad options\n";
die "usage: $0 [--tol=T] [--floor=F] file1 file2\n" unless @ARGV == 2;
sub read_xy {
  my $fn = shift;
  open(my $f, "<", $fn) or die "cannot read $fn: $!\n";
  my @d;
  while (<$f>) { next if /^\s*(#|$)/; push @d, [split]; }
  close $f;
  return \@d;
}
my ($a, $b) = map { read_xy($_) } @ARGV;
die "different numbers of points: " . scalar(@$a) . " and " . scalar(@$b) . "\n" unless @$a == @$b;
my ($worst, $bad) = (0, 0);
for my $i (0 .. $#$a) {
  die "different mesh at point $i\n" unless abs($a->[$i][0] - $b->[$i][0]) <= 1e-12 * abs($a->[$i][0]);
  my $dev = abs($a->[$i][1] - $b->[$i][1]);
  my $allowed = $tol * abs($a->[$i][1]) + $floor;
  $worst = $dev / $allowed if $dev / $allowed > $worst;
  if ($dev > $allowed) { print "point $i x=$a->[$i][0]: $a->[$i][1] vs $b->[$i][1]\n" if $bad < 10; $bad++; }
}
printf "%d points, largest deviation %.3g of the allowed\n", scalar(@$a), $worst;
exit($bad ? 1 : 0);
