# Get Ancestry proportions for each hybrid based on admixfrog

```
for f in *concatforperms.txt; do
    ./Calculate_ancestry_proportions.pl "$f"
done > ancestry_summary.txt
```

Perl script:
```
#!/usr/bin/perl

use strict;
use warnings;

my $file = shift or die "Usage: $0 file.csv\n";

open(my $fh, "<", $file) or die "Cannot open $file: $!\n";

my %sum;
my @ancestries;
my %idx;
my $n = 0;

while (<$fh>) {

    chomp;
    next unless /\S/;

    # Process header lines (including repeated headers)
    if (/^chrom,/) {

        my @header = split /,/;

        %idx = ();
        @idx{@header} = (0 .. $#header);

        @ancestries = ();

        my $after_nsnps = 0;

        for my $col (@header) {

            if ($col eq "n_snps") {
                $after_nsnps = 1;
                next;
            }

            next unless $after_nsnps;

            if ($col =~ /^[A-Z]{3}$/) {
                push @ancestries, $col;
            }
        }

        next;
    }

    my @f = split /,/;

    my %row;

    # Pure ancestry contributions
    for my $anc (@ancestries) {
        $row{$anc} = $f[$idx{$anc}];
    }

    # Split pairwise ancestry probabilities evenly
    for my $i (0 .. $#ancestries - 1) {

        for my $j ($i + 1 .. $#ancestries) {

            my $a = $ancestries[$i];
            my $b = $ancestries[$j];

            my $pair;

            if (exists $idx{"$a$b"}) {
                $pair = "$a$b";
            }
            elsif (exists $idx{"$b$a"}) {
                $pair = "$b$a";
            }
            else {
                next;
            }

            my $val = $f[$idx{$pair}];

            $row{$a} += 0.5 * $val;
            $row{$b} += 0.5 * $val;
        }
    }

    for my $anc (@ancestries) {
        $sum{$anc} += $row{$anc};
    }

    $n++;
}

die "No data rows found\n" if $n == 0;

print $file;

my $total = 0;

for my $anc (@ancestries) {

    my $avg = $sum{$anc} / $n;

    $total += $avg;

    printf " %s:%.6f", $anc, $avg;
}

printf " TOTAL:%.6f\n", $total;

close($fh);
```
