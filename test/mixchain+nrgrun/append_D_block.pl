#!/usr/bin/env perl
# Append the D block of adaptive broadening to a data file that a template's instantiate wrote: the chain at
# z' = z + dz from the staged coefficient files in DIRECTORY (mixchain's dz/, after stage_<SYM>.pl), in the layout of
# the template's z block, which is the count line and the table of each xi set, then the same for each zeta set.
#
#   append_D_block.pl SETS [DIRECTORY [DATA]]      (defaults: dz, data)
#
# SETS is the number of coefficient sets: 1 for QS, 4 for U1 with pol2x2. The shift is read from DIRECTORY/dz.dat,
# and the count line is taken from the z block of DATA, so that the two blocks have the same shape by construction;
# every table must hold count+1 values.
use strict;
use warnings;

my $sets      = shift // die "usage: append_D_block.pl SETS [DIRECTORY [DATA]]\n";
my $directory = shift // "dz";
my $data      = shift // "data";
die "append_D_block: SETS must be a positive integer\n" unless $sets =~ /^[1-9]\d*$/;

sub read_lines {
  my $filename = shift;
  open(my $in, "<", $filename) or die "append_D_block: cannot read $filename: $!\n";
  my @lines = grep { !/^\s*(#|$)/ } <$in>;
  close $in;
  chomp @lines;
  return @lines;
}

my @dz = read_lines("$directory/dz.dat");
die "append_D_block: $directory/dz.dat must hold one number\n" unless @dz == 1;

# The count line of the z block: the line after the one holding only "z".
my @lines = read_lines($data);
my ($count, $seen) = (undef, 0);
for my $i (0 .. $#lines - 1) {
  next unless $lines[$i] =~ /^z\s*$/;
  die "append_D_block: $data has more than one z block\n" if $seen++;
  $count = $lines[$i + 1];
}
die "append_D_block: $data has no z block\n" unless defined $count;
die "append_D_block: $data has a D block already\n" if grep { /^D / } @lines;
$count =~ s/^\s+|\s+$//g;

open(my $out, ">>", $data) or die "append_D_block: cannot append to $data: $!\n";
print $out "D $dz[0]\n";
for my $kind ("xi", "zeta") {
  for my $set (1 .. $sets) {
    my @values = read_lines("$directory/$kind$set.dat");
    die "append_D_block: $directory/$kind$set.dat holds " . scalar(@values) . " values, the z block " . ($count + 1) . "\n"
      unless @values == $count + 1;
    print $out "$count\n";
    print $out "$_\n" for @values;
  }
}
close $out or die "append_D_block: error writing $data: $!\n";
