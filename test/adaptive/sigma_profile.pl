#!/usr/bin/env perl
# Profile of the broadening widths in an .abin file (adaptive broadening): per decade of |omega|, the total weight,
# the weight-averaged sigma/ln(Lambda) and the weighted 10%, 50% and 90% quantiles of sigma/ln(Lambda).
#
#   sigma_profile.pl --lambda=L [--columns=3|4] [--check [--step=0.01]] file.abin
#
# .abin: headerless native doubles, records (omega, sigma, Re w[, Im w]); 3 columns unless reim=true.
#
# --check applies the conditions of test/adaptive/adaptive_sigma (a SIAM with a flat band, T = 0+) and exits with 1
# if one fails: the total weight is 1 within 1e-6; in every decade from 1e-9 to 1e-4 (the strong-coupling fixed
# point) <sigma> is within 1% of the sigma grid point nearest ln(Lambda); <sigma>/ln(Lambda) < 0.85 in [1, 10) (the
# band edge) and < 0.97 in [0.1, 1) (the Hubbard peaks).
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(floor);

my ($lambda, $columns, $check, $step) = (undef, 3, 0, 0.01);
GetOptions("lambda=f" => \$lambda, "columns=i" => \$columns, "check" => \$check, "step=f" => \$step) or die "bad options\n";
die "usage: $0 --lambda=L [--columns=3|4] file.abin\n" unless defined $lambda && @ARGV == 1;
my $fn = $ARGV[0];
open(my $f, "<:raw", $fn) or die "cannot read $fn: $!\n";
local $/;
my $raw = <$f>;
close $f;
my @v = unpack("d*", $raw);
die "$fn: size not a multiple of $columns doubles\n" if @v % $columns;
my $lnl = log($lambda);

my (%cells, $total, $nonfinite);
for (my $r = 0; $r < @v; $r += $columns) {
  my ($e, $s, $w) = @v[$r .. $r + 2];
  if ($e != $e || $s != $s || $w != $w || abs($e) == 9**9**9 || abs($w) == 9**9**9) { $nonfinite++; next; }
  my $decade = floor(log(abs($e)) / log(10));
  push @{$cells{$decade}}, [$s / $lnl, $w];
  $total += $w;
}
die "$fn: $nonfinite non-finite records\n" if $nonfinite;
printf "%s: %d cells, total weight %.6f\n", $fn, @v / $columns, $total;
printf "%-18s %12s %9s %7s %7s %7s\n", "|omega| decade", "weight", "<s>/lnL", "q10", "q50", "q90";
my %mean; # decade -> <sigma>/ln(Lambda)
for my $d (sort { $b <=> $a } keys %cells) {
  my @c = sort { $a->[0] <=> $b->[0] } @{$cells{$d}};
  my ($sw, $sws, $sabs) = (0, 0, 0);
  for (@c) { $sw += $_->[1]; $sws += $_->[1] * $_->[0]; $sabs += abs($_->[1]); }
  next if $sabs == 0;
  my $quantile = sub {
    my $q = shift;
    my $acc = 0;
    for (@c) { $acc += abs($_->[1]); return $_->[0] if $acc >= $q * $sabs; }
    return $c[-1][0];
  };
  $mean{$d} = $sw != 0 ? $sws / $sw : 0;
  printf "[1e%+03d, 1e%+03d) %12.4e %9.4f %7.3f %7.3f %7.3f\n", $d, $d + 1, $sw,
    $mean{$d}, $quantile->(0.1), $quantile->(0.5), $quantile->(0.9);
}

exit 0 unless $check;
my @failures;
push @failures, sprintf("total weight %.9f, not 1 within 1e-6", $total) if abs($total - 1) > 1e-6;
my $fixed_point = $step * POSIX::floor($lnl / $step + 0.5) / $lnl; # the grid point nearest ln(Lambda)
for my $d (-9 .. -5) {
  if (!exists $mean{$d}) { push @failures, "no weight in the decade 1e$d"; next; }
  push @failures, sprintf("decade 1e%d: <s>/lnL = %.4f, not within 1%% of %.4f", $d, $mean{$d}, $fixed_point)
    if abs($mean{$d} - $fixed_point) > 0.01 * $fixed_point;
}
push @failures, sprintf("band edge [1, 10): <s>/lnL = %.4f, not below 0.85", $mean{0} // -1) unless defined $mean{0} && $mean{0} < 0.85;
push @failures, sprintf("Hubbard peaks [0.1, 1): <s>/lnL = %.4f, not below 0.97", $mean{-1} // -1) unless defined $mean{-1} && $mean{-1} < 0.97;
if (@failures) { print "FAIL $_\n" for @failures; exit 1; }
printf "OK (fixed point %.4f, band edge %.4f, Hubbard peaks %.4f)\n", $fixed_point, $mean{0}, $mean{-1};
