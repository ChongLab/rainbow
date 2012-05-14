#!/usr/bin/perl -w
#

use strict;
use warnings;

my $file = shift or die "Usage: $0 <rbasmed file>\n";

my $len = 0;
my $name = "";
my $seq = "";
open IN, $file or die $!;
while (<IN>) {
	if (/^E/) {
		my @e = split /\s+/, $_;
		if ($len != 0) {
			print ">$name"."_L"."$len\n";
			print $seq, "\n";
		}
		$name = $e[0].$e[1];
		$len = 0;
		$seq = "";
	} elsif (/^S/) {
		my @e = split /\s+/, $_;
		if ($len < length($e[1])){
			$seq = $e[1];
			$len = length $e[1];
		}
	}
}
close IN;

print ">$name"."_L"."$len\n";
print $seq, "\n";
