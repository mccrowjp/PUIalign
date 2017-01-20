#!/usr/bin/env perl
#
# PUIalign - Phylogenetically Unambiguous Indel Alignment
#
# Created by John P. McCrow - 9/5/2007
#
# Citation:
#    McCrow JP, Alignment of phylogenetically unambiguous indels in Shewanella,
#    Journal of Computational Biology 16 (11), 1517-1528 (2009)
#

use strict;
use progress;

my $indir = shift;
my $infile = shift;

my $matchlen = 1;
my $alignlen = 20;
my $minindellen = 1;
my $minspecies = 4;

my @filelist;
my %protsp;

unless($indir && $infile) {
    die "Usage: $0 [Alignments Directory] [Prot Species Table]\n";
}

foreach my $file (glob($indir."/*")) {
    push(@filelist, $file);
}

open(IN, $infile) or die "Unable to open file $infile\n";

while(<IN>) {
    chomp;
    unless(/^\#/) {
        my ($prot, $sp) = split(/\t/);
        $protsp{$prot}=$sp;
    }
}
close(IN);

my $numdone = 0;
my $starttime = gettimeofday();

foreach my $file (@filelist) {
    drawprogress($numdone, scalar(@filelist), $starttime);
    $numdone++;

    open(IN, $file) or die "Unable to open file $file\n";

    my %refstr = ();
    my %refsp = ();
    my %seqpos = ();
    my %refpos = ();
    my %delpos = ();
    my %allpos = ();

    while(<IN>) {
        chomp;
        if(/^[^\s]+\s+[^\s]+$/) {
            my ($prot, $str) = split(/\s+/);
            my $protid = $prot;
            
            if($prot =~ /\|.+\|.+\|/) {
                $protid = (split(/\|/, $prot))[3];
            }
            $refstr{$protid} .= $str;
        }
    }
    close(IN);

    #Remove possible paralogs
    foreach my $r (keys %refstr) {
        $refsp{$protsp{$r}}++;
    }
    foreach my $r (keys %refstr) {
        if($refsp{$protsp{$r}} > 1) {
            delete($refstr{$r});
        }
    }

    if(scalar(keys %refstr) >= $minspecies) {

        #Map multiple alignment positions to sequence positions for all sequences
        foreach my $r (keys %refstr) {
            my $j=0;
            for(my $i=0; $i<length($refstr{$r}); $i++) {
                $seqpos{$r.",".$i} = $j;
                if(substr($refstr{$r},$i,1) ne "-") {
                    $j++;
                } 
            }
        }

        #Find deletions in multiple alignments
        foreach my $r (keys %refstr) {
            while($refstr{$r} =~ /\w{$matchlen}(\-{$minindellen,})\w{$matchlen}/gi) {
                my $curdelpos = ($-[1]).",".($+[1]-1);
                #If $alignlen length sequence on either side than include this indel
                if($seqpos{$r.",".($-[1])} > $alignlen 
                   && $seqpos{$r.",".($+[1])} < $seqpos{$r.",".(length($refstr{$r})-1)} - $alignlen) {
                    $refpos{$r.",".$curdelpos} = 1;
                    $allpos{$curdelpos}++;
                }
            }
        }

        if(scalar(keys %allpos) > 0) {

            #Print sequences and potential indels
            print ">$file\n";
            foreach my $r (sort keys %refstr) {
                my $s = $refstr{$r};
                $s =~ s/-//g;
                print join("\t", ($protsp{$r}, $r, $s))."\n";
            }
            
            my $lasta2 = -1;
            foreach my $a (sort {$a<=>$b} keys %allpos) {
                my ($a1,$a2) = split(/,/, $a);
                if($a1 > $lasta2) {
                    my $a3 = int(($a1 + $a2) / 2);
                    my @poslist = ();
                    foreach my $r (sort keys %refstr) {
                        push(@poslist, $seqpos{$r.",".$a3});
                    }
                    print "#".($a2-$a1+1).":".join(",", @poslist)."\n";
                    $lasta2 = $a2;
                }
            }
            print "\n";
        }
    }
}

eraseprogress();
