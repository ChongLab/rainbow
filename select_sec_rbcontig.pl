#!/usr/bin/perl -w
#

use strict;
use warnings;

my $file = shift or die "Usage: $0 <rbasmed file>\n";

my $len = 0;
my $name = "";
my $seq = "";
my $seclen = 0;
my $secseq = "";
open IN, $file or die $!;
while (<IN>) {
	if (/^E/) {
		my @e = split /\s+/, $_;
		if ($len) {
			print ">$name"."_L"."$len\n";
			print $seq, "\n";
		}
		if ($seclen) {
			print ">$name"."_L"."$seclen\n";
			print $secseq, "\n";
		}
		$name = $e[0].$e[1];
		$len = $seclen = 0;
		$seq = $secseq = "";
	} elsif (/^S/) {
		my @e = split /\s+/, $_;
		if ($len < length($e[1])){
			$secseq = $seq;
			$seq = $e[1];
			$seclen = length $secseq;
			$len = length $e[1];
		} elsif (length($e[1]) > $seclen) {
			$secseq = $e[1];
			$seclen = length $secseq;
		}
	}
}
close IN;

print ">$name"."_L"."$len\n";
print $seq, "\n";
if ($seclen) {
	print ">$name"."_L"."$seclen\n";
	print $secseq, "\n";
}
