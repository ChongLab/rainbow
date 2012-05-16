#!/usr/bin/perl -w
#
use strict;

my $merged = shift or &usage();
my $div = shift or &usage();
my $asm = shift or &usage();


my %mhash = ();

open IN, $merged or die $!;
while (<IN>) {
	chomp;
	my @e = split /\s+/, $_;
	if (@e > 1) {
		foreach (@e) {
			$mhash{$_} = $e[0];
		}
	}
}
close IN;

open OUT, ">tmp_to_reasm.rbout" or die $!;
open IN, $div or die $!;
while (<IN>) {
	chomp;
	my @e = split /\s+/, $_;
	if ($mhash{$e[1]}) {
		$e[1] = $mhash{$e[1]};
		print OUT join("\t", @e), "\n";
	}
}
close IN;

close OUT;

my $len = 0;
my $name = "";
my $seq = "";
my $pre = "";
open OUT, ">final_asm.fa" or die $!;
open IN, $asm or die $!;
my @e = ();
while (<IN>) {
	if (/^E/) {
		@e = split /\s+/, $_;
		if ($len && !$mhash{$pre}) {
			print OUT ">$name"."_L"."$len\n";
			print OUT $seq, "\n";
		}
		$name = $e[0].$e[1];
		$len = 0;
		$seq = "";
		$pre = $e[1];
	} elsif (/^S/) {
		@e = split /\s+/, $_;
		if ($len < length $e[1]) {
			$seq = $e[1];
			$len = length $e[1];
		}
	}
}
close IN;
if ($len && !$mhash{$pre}) {
	print OUT ">$name"."_L"."$len\n";
	print OUT $seq, "\n";
}
close OUT;

`sort -k2,2n tmp_to_reasm.rbout| tee tmp_to_reasm.rbout.srt |rbasm -i - |  select_best_rbcontig.pl - | sed 's/^>E/>R/' >> final_asm.fa`;
`rm tmp_to_reasm.rbout*`;

1;

sub usage {
	die(qq/
Usage: rerun_rbsam.pl <merged out> <divided out> <rbasm out>
       OUTPUT: final_asm.fa
       NOTE: program 'rbasm' and 'select_best_rbcontig.pl' must be
       copied to the 'bin' or set in the correct 'PATH'
\n/);
}
