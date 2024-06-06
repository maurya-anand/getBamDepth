#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

my ($input_bed, $input_bam, $input_depth, $thresholds_str);

GetOptions(
    'bed=s'       => \$input_bed,
    'bam=s'       => \$input_bam,
    'depth=s'     => \$input_depth,
    'thresholds=s' => \$thresholds_str,
) or die "Usage: $0 --bed BED_FILE [--bam BAM_FILE | --depth DEPTH_FILE] [--thresholds THRESHOLDS]\n";

die "Bed file is required\n" unless $input_bed;
die "Bed file does not exist: $input_bed\n" unless -e $input_bed;

die "Either bam or depth file is required\n" unless $input_bam || $input_depth;
die "Bam file does not exist: $input_bam\n" if $input_bam && !-e $input_bam;
die "Depth file does not exist: $input_depth\n" if $input_depth && !-e $input_depth;

my @thresholds = $thresholds_str ? split /,/, $thresholds_str : (10, 20, 30, 40, 50, 60, 70, 80, 90, 100);

my $sample;
if ($input_bam) {
    $sample = fileparse($input_bam, qr/\.[^.]*/);
    $input_depth = "$sample.depth";
    `samtools depth $input_bam > $input_depth`;
} elsif ($input_depth) {
    $sample = fileparse($input_depth, qr/\.[^.]*/);
}

my @header = ("sample", "chrom", "start", "end", "region", "avg_depth", map { $_ . "x" } @thresholds);
my @results;

open (BED, $input_bed) or die "Can't read the bed file: $input_bed";
while (my $rec = <BED>) {
    chomp $rec;
    my ($chr, $start, $end, $info, $temp, $strand) = split(/\t/, $rec);
    $start += 1;
    my @depths = get_depth($chr, $start, $end, $input_depth, @thresholds);
    push @results, ["$sample", $chr, $start, $end, $info // "unknown", @depths];
}
close BED;

print join("\t", @header) . "\n";
foreach my $result (@results) {
    print join("\t", @$result) . "\n";
}

system("rm $input_depth");

sub get_depth {
    my ($chr, $start, $end, $depth, @thresholds) = @_;
    open (DEP, $depth) or die "Can't open $depth: $!";
    
    my ($sum, $count) = (0, 0);
    my @counters = (0) x scalar @thresholds;

    while (my $line = <DEP>) {
        chomp $line;
        my @fields = split /\t/, $line;
        if ($fields[0] eq $chr && $fields[1] >= $start && $fields[1] <= $end) {
            my $depth_val = $fields[2];
            $sum += $depth_val;
            $count++;
            for my $i (0..$#thresholds) {
                if ($depth_val >= $thresholds[$i]) {
                    $counters[$i]++;
                }
            }
        }
    }
    close DEP;
    
    my $avg_depth = $count > 0 ? sprintf("%.2f", $sum / $count) : "0";
    return ($avg_depth, @counters);
}

