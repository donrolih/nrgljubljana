#!/usr/bin/env perl
# Compare dE/dz of the middle run with central finite differences of the energies of three runs at z-dz, z, z+dz.
#
#   compare_hf.pl --lambda=L --dz=DZ --tol=TOL minus/dEdz.nrg middle/dEdz.nrg plus/dEdz.nrg
#
# The files are written with adaptive_broadening_dumpenergies: per iteration and subspace the energies E (in units of
# the shell, which depends on z through SCALE ~ Lambda^-z) and dE/dz (physical derivative in units of the shell,
# relative to the ground state). In units of the shell of the middle run, the physical derivative of E_i - E_GS is
#
#   -ln(Lambda) (E_i - E_GS) + [(E_i - E_GS)(z+dz) - (E_i - E_GS)(z-dz)] / (2 dz),
#
# the first term from the shell unit, the second from the levels themselves. Each state must agree within
# TOL * max(1, |expected|). All three runs must have the same subspaces and numbers of states (no truncation).
#
# Near an avoided crossing within [z-dz, z+dz] the difference quotient of a single level is meaningless, while
# Hellmann-Feynman still holds. States of a subspace whose neighbouring levels come closer than GAP in any of the
# three runs are therefore grouped, and a group is compared through its sum: the trace of dH/dz over the group is the
# derivative of the summed energies, which is smooth through the crossing.
use strict;
use warnings;
use Getopt::Long;

my ($lambda, $dz, $tol, $gap) = (undef, undef, undef, 1e-3);
GetOptions("lambda=f" => \$lambda, "dz=f" => \$dz, "tol=f" => \$tol, "gap=f" => \$gap) or die "bad options\n";
die "usage: $0 --lambda=L --dz=DZ --tol=TOL [--gap=GAP] minus middle plus\n" unless defined $lambda && defined $dz && defined $tol && @ARGV == 3;

# Returns { iteration => { subspace => { E => [...], dEdz => [...] } } }
sub read_dump {
  my $fn = shift;
  open(my $f, "<", $fn) or die "cannot read $fn: $!\n";
  my (%d, $n, $sub);
  while (my $line = <$f>) {
    if ($line =~ /^===== Iteration number: (-?\d+)/) {
      $n = $1;
      die "$fn: iteration $n appears twice\n" if exists $d{$n};
      $d{$n} = {};
      next;
    }
    if ($line =~ /^Subspace: (.*\S)\s*$/) { $sub = $1; next; }
    if ($line =~ /^(E|dEdz)= (.*)$/) {
      die "$fn: values before a subspace\n" unless defined $n && defined $sub;
      $d{$n}{$sub}{$1} = [split ' ', $2];
    }
  }
  close $f;
  die "$fn: no data\n" unless %d;
  return \%d;
}

# E_i - E_GS for one iteration, GS = the lowest energy over all subspaces.
sub relative {
  my $it = shift;
  my $gs;
  for my $s (values %$it) { for my $e (@{$s->{E}}) { $gs = $e if !defined $gs || $e < $gs; } }
  my %r;
  for my $sub (keys %$it) { $r{$sub} = [map { $_ - $gs } @{$it->{$sub}{E}}]; }
  return \%r;
}

my ($m, $o, $p) = map { read_dump($_) } @ARGV;
my $lnl = log($lambda);
my ($count, $groups, $worst, $where) = (0, 0, 0.0, "");
my $failed = 0;
for my $n (sort { $a <=> $b } keys %$o) {
  die "iteration $n missing in a neighbouring run\n" unless $m->{$n} && $p->{$n};
  my ($rm, $ro, $rp) = (relative($m->{$n}), relative($o->{$n}), relative($p->{$n}));
  for my $sub (sort keys %{$o->{$n}}) {
    my $d = $o->{$n}{$sub}{dEdz} or die "iteration $n subspace $sub: no dEdz\n";
    die "iteration $n subspace $sub: different number of states in the three runs\n"
      unless $rm->{$sub} && $rp->{$sub} && @{$rm->{$sub}} == @$d && @{$rp->{$sub}} == @$d && @{$ro->{$sub}} == @$d;
    # Groups of consecutive states (sorted by energy) closer than $gap in any run.
    my @group = (0);
    my @groups_here;
    for my $i (1 .. $#$d) {
      my $close = 0;
      for my $r ($rm, $ro, $rp) { $close = 1 if abs($r->{$sub}[$i] - $r->{$sub}[$i-1]) < $gap; }
      if ($close) { push @group, $i; } else { push @groups_here, [@group]; @group = ($i); }
    }
    push @groups_here, [@group] if @$d;
    for my $g (@groups_here) {
      my ($got, $expected) = (0.0, 0.0);
      for my $i (@$g) {
        $got += $d->[$i];
        $expected += -$lnl * $ro->{$sub}[$i] + ($rp->{$sub}[$i] - $rm->{$sub}[$i]) / (2 * $dz);
      }
      my $dev = abs($got - $expected) / (abs($expected) > 1 ? abs($expected) : 1);
      my $label = @$g == 1 ? "state $g->[0]" : "states $g->[0]-$g->[-1] (summed)";
      $count += @$g;
      $groups++ if @$g > 1;
      if ($dev > $worst) { $worst = $dev; $where = "iteration $n subspace $sub $label: dEdz=$got expected=$expected"; }
      if ($dev > $tol) {
        print "FAIL iteration $n subspace $sub $label: dEdz=$got expected=$expected\n" if $failed < 20;
        $failed++;
      }
    }
  }
}
printf "%d states compared (%d groups of nearly degenerate states by their sums), largest deviation %.3e (%s)\n",
  $count, $groups, $worst, $where;
if ($failed) { print "$failed states outside the tolerance $tol\n"; exit 1; }
print "OK\n";
