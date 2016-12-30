package Gene_util;
use strict;
use warnings;
use POSIX;

use Exporter qw(import);

our @EXPORT_OK = qw(system2 read_header);

my $bin = "gene_bin.tar.gz";
my $delim = "\x01";

sub get_gene_bin {
    return $bin;
}

sub get_field_delim {
    return $delim;
}

my @std_ref= (
    "1" => 249250621,
    "2" => 243199373,
    "3" => 198022430,
    "4" => 191154276,
    "5" => 180915260,
    "6" => 171115067,
    "7" => 159138663,
    "8" => 146364022,
    "9" => 141213431,
    "10" => 135534747,
    "11" => 135006516,
    "12" => 133851895,
    "13" => 115169878,
    "14" => 107349540,
    "15" => 102531392,
    "16" => 90354753,
    "17" => 81195210,
    "18" => 78077248,
    "19" => 59128983,
    "20" => 63025520,
    "21" => 48129895,
    "22" => 51304566,
    "X" => 155270560,
    "Y" => 59373566,
    "MT" => 16569);

sub get_std_ref {
    return @std_ref;
}

my @seq2rname = (1 .. 22, "X", "Y", "MT");

sub get_chr_seq {
    my %chr_seq;
    for (my $i = 0; $i <= $#seq2rname; $i++) {
        $chr_seq{$seq2rname[$i]} = $i + 1;
    }
    return %chr_seq;
}

sub check_eof {
    if (eof(STDIN)) {
        print STDERR "no input, exit now.\n";
        exit;
    }
}

sub system2 {
    my $cmd = shift;
    my $t = time;
    print STDERR strftime("%F %T", localtime), " >> start command: $cmd\n";
    0 == system($cmd) or die "command failed: $cmd\n";
    $t = time - $t;
    print STDERR strftime("%F %T", localtime), " << command finished in $t seconds: $cmd\n";
}

sub read_header {
    my $header_file = shift;
    open my $fh, '<', $header_file;
    my $header = <$fh>;
    close $fh;
    return $header;
}

1;
