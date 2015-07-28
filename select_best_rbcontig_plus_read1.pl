#!/usr/bin/perl -w
#

use strict;
use warnings;

my $file = shift or die "Usage: $0 <rbasmed file> <rbdiv output file>\n";
my $div = shift or die "Usage: $0 <rbasmed file> <rbdiv output file>\n";

my $len = 0;
my $name = "";
my $seq = "";
my $reads = "";

open DIV, $div or die $!;
open IN, $file or die $!;
my %cls = ();
$_ = <DIV>;
my $dseq = $_;
my @de = split /\s+/, $dseq;
while (<IN>) {
	if (/^E/) {
		my @e = split /\s+/, $_;
		if ($len != 0) {
			if (&isoverlap($reads)) {
				print ">$name"."_L"."$len\n";
				print $seq, "\n"
			} else {
				print ">$name"."_L"."$len\n";
				print &gen_mock_ctg($seq, $de[2]), "\n";
			}
			while (<DIV>) {
				my @e2 = split /\s+/, $_;
				if ($e2[1] eq $e[1]) { # thanks to Jonathan Puritz for pointing this
					$dseq = $_;
					@de = @e2;
					last;
				}
			}
		}
		$name = $e[0].$e[1];
		$len = 0;
		$seq = "";
		$reads = "";
	} elsif (/^S/) {
		my @e = split /\s+/, $_;
		if ($len < length($e[1])) {
			$seq = $e[1];
			$len = length $e[1];
			<IN>;
			$_ = <IN>;
			$reads = $_;
		}
	} 
}
close IN;
close DIV;
if (&isoverlap($reads)) {
	print ">$name"."_L"."$len\n";
	print $seq, "\n"
} else {
	print ">$name"."_L"."$len\n";
	print &gen_mock_ctg($seq, $de[2]), "\n";
}

1;

sub isoverlap {
	my $reads = shift;
	my @ids = split /\s+/, $reads;
	my $ret = 0;
	foreach my $id (@ids) {
		next if $id =~ /R/;
		my @rec = split /:/, $id;
		if ($rec[2] == 0) {
			$ret = 1;
			return $ret;
		}
	}
	return $ret;
}

sub gen_mock_ctg {
	my @seqs = @_;
	$seqs[1] =~ tr/ACGTacgt/TGCAtgca/;
	$seqs[1] = reverse $seqs[1];
	return $seqs[0].("N"x10).$seqs[1];
}
