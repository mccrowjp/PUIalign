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

my ($indir, $inall, $inindels, $inps, $intree) = @ARGV;

my $viewsize = 40;
my $excluderadius = 20;
my $maxindelgroups = 3;
my $minindellen = 1;
my $mintaxa = 4;
my $alpha = 2;
my $beta = 10;

my %spname;
my %spnum;
my %protsp;
my %allspecies;
my %allspecies_ps;
my %allspecies_i;
my %allspecies_t;
my %groupingname;
my %grouping;
my %numgroups;
my %allgenes;
my %pboundl;
my %pboundu;
my %minbeta;
my %minscore;
my %goodgenes;
my %goodgene_ph;
my %numtaxa;
my %indelgroupprots;
my %indelprotgroup;
my %indelgapsize;
my %indellens;
my %indelprotpos;
my %indelprotlist;
my %geneindelpos;
my %exclude;
my %alphacutcount;
my %branchcount;
my %branchcountunique;
my %branchcount_good;
my %branchcountunique_good;
my %phstr;
my %splitstr;
my %splitstrincl;

my $numindels;
my $numgoodindels;
my $numphchrs;
my $maxsplen;
my $numintree;
my $numouttree;
my $numintree_ph;
my $numouttree_ph;


unless($indir && $inall && $inindels && $inps) {
    die "Usage: $0 [Alignments Directory] [All Indels] [Indel Group Output] [Prot Species Table] ([Tree Groups])\n";
}

print "#Alignments Dir: $indir\n";
print "#All Indels: $inall\n";
print "#Indel Groups: $inindels\n";
print "#Prot-Species: $inps\n";
if($intree) {
    print "#Tree Groups: $intree\n";
}
print "#Viewsize = $viewsize\n";
print "#Max Indel Groups = $maxindelgroups\n";
print "#Min Gap Length = $minindellen\n";
print "#Alpha = $alpha\n";
print "#Beta = $beta\n";
print "\n";

open(INS, $inps) or die "Unable to open file $inps\n";
while(<INS>) {
    chomp;
    if(/^\#(\d+)\s+([^\s]+)\s*$/) {
        $spname{$1}=$2;
        $spnum{$2}=$1;	
    } else {
        my ($prot, $sp) = split(/\t/);
        $protsp{$prot}=$sp;
        $allspecies_ps{$sp}=1;
        if(length($spname{$sp}) > $maxsplen) {
            $maxsplen = length($spname{$sp});
        }
    }
}
close(INS);

my $numgroupings = 0;
if($intree) {
    open(INT, $intree) or die "Unable to open file $intree\n";
    while(<INT>) {
        chomp;
        if(length($_)>1) {
            print "$numgroupings: $_\n";
            $groupingname{$numgroupings} = $_;
            my @groupsp = split(/,/);
            foreach my $s (@groupsp) {
                if($allspecies_ps{$spnum{$s}}) {
                    $grouping{"".$numgroupings.",".$spnum{$s}} = 1;
                    $allspecies_t{$spnum{$s}} = 1;
                }
            }
            $numgroupings++;
        }
    }
    close(INT);
    print "\n";
}

open(INI, $inindels) or die "Unable to open file $inindels\n";
my $curindelname;
my $curgene;

while(<INI>) {
    chomp;
    unless(/^\#/) {
        if(/^>(.+)$/) {
            $curindelname = $1;
            $numgroups{$curindelname}=0;
            ($curgene) = split(/,/, $curindelname);
            $allgenes{$curgene}=1;
            $numindels++;
        } elsif(/^(\d+),(\d+),(\d+),(\d+)$/) {
            $pboundl{$curindelname} = $1;
            $pboundu{$curindelname} = $2;
            $minbeta{$curindelname} = $3;
            $minscore{$curindelname} = $4;
            $goodgenes{$curgene}=1;
        } elsif(/^\[([^\]]+)\]:(.+)$/) {
            $numgroups{$curindelname}++;
            my @protids = split(/,/, $2);
            $numtaxa{$curindelname} += scalar(@protids);
            $indelgroupprots{$curindelname.":".$numgroups{$curindelname}} = join(",",@protids);
            foreach my $p (@protids) {
                $indelprotgroup{$curindelname.":".$p} = $numgroups{$curindelname};
                $allspecies_i{$protsp{$p}}=1;
            }

            foreach my $l (split(/,/,$1)) {
                if(!defined($indelgapsize{$curindelname}) || abs($l) > $indelgapsize{$curindelname}) {
                    $indelgapsize{$curindelname} = abs($l);
                }
                $indellens{abs($l)}++;
            }
        }
    }
}
close(INI);



#Only show species that have data in indels, protein_species, and if applicable, tree groups
foreach my $s (%allspecies_i) {
    if($allspecies_i{$s} && $allspecies_ps{$s} && (!$intree || $allspecies_t{$s})) {
        $allspecies{$s} = 1;
    }
}

open(INA, $inall) or die "Unable to open file $inall\n";
my $curgene;
my $indelnum;

while(<INA>) {
    chomp;
    if(/^>(.+)$/) {
        $curgene = $1;
        $indelnum = 1;
    } elsif(/^\#\d+:(.+)$/) {
        my @spposlist = split(/,/, $1);
        for(my $i=0; $i<scalar(@spposlist); $i++) {
            $indelprotpos{$curgene.",".$indelnum.":".@{$indelprotlist{$curgene}}[$i]} = $spposlist[$i];
        }
        $indelnum++;
    } else {
        my ($sp, $prot, $str) = split(/\t/);
        push(@{$indelprotlist{$curgene}}, $prot);	
    }
}
close(INA);

my $starttime = gettimeofday();
my $totindels = scalar(keys %numgroups);
my $totdone = 0;

foreach my $indel (sort {($pboundu{$b}-$pboundl{$b})<=>($pboundu{$a}-$pboundl{$a})} keys %numgroups) {
    drawprogress($totdone, $totindels, $starttime);
    $totdone++;

    my %branches = ();
    
    if($pboundu{$indel} && $numgroups{$indel} >= 2 && $numgroups{$indel} <= $maxindelgroups) {
        my ($gene,$inum) = split(/,/, $indel);
        my $file = $gene;

        if($intree) {
            #Fit indel to branches by each group
            for(my $g1=1; $g1<=$numgroups{$indel}; $g1++) {
                my %group1 = ();
                my %group2 = ();

                #Group1 is all species in one group of this indel
                foreach my $p (split(/,/, $indelgroupprots{$indel.":".$g1})) {
                    $group1{$protsp{$p}} = 1;
                }
                #Group2 is all the other species in other groups of this indel
                for(my $g2=1; $g2<=$numgroups{$indel}; $g2++) {
                    if($g2 != $g1) {
                        foreach my $p (split(/,/, $indelgroupprots{$indel.":".$g2})) {
                            $group2{$protsp{$p}} = 1;
                        }
                    }
                }

                my @found = ();
                for(my $i=0; $i<$numgroupings; $i++) {
                    my $a = "";
                    my $b = "";
                    foreach my $s (keys %allspecies) {
                        if($group1{$s}) {
                            $a .= (0+$grouping{"".$i.",".$s});
                        } elsif($group2{$s}) {
                            $b .= (0+$grouping{"".$i.",".$s});
                        }
                    }
                    
                    #$allin = ($a.$b) =~ /^1+$/;
                    my $opposing = (($a =~ /0/ && $b =~ /1/) || ($a =~ /1/ && $b =~ /0/));
                    my $aisgood = !($a =~ /01/ || $a =~ /10/);
                    my $bisgood = !($b =~ /01/ || $b =~ /10/);

                    #print "$i $a:$b\n";

                    if($opposing && ($aisgood && $bisgood)) {
                        push(@found, $i);
                    }
                }

                if(scalar(@found) > 0) {
                    foreach my $b (@found) {
                        $branches{$b} = 1;
                    }
                }
            }
        }

        open(IN, $indir."/".$file) or die "Unable to open file $indir/$file\n";
        my %refstr = ();

        while(<IN>) {
            chomp;
            if(/^[^\s]+[\s]+[^\s]+$/) {
                my ($prot, $str) = split(/\s+/);
                my $protid = $prot;
                if($prot =~ /\|.+\|.+\|/) {
                    $protid = (split(/\|/, $prot))[3];
                }
                $refstr{$protid} .= $str;
            }
        }
        close(IN);

        #Find left position in multiple alignment
        my $r = @{$indelprotlist{$gene}}[0];
        my $j = $indelprotpos{$indel.":".@{$indelprotlist{$gene}}[0]};
        my $i = 0;
            
        while($j>0 && $i<length($refstr{$r})) {
            if(substr($refstr{$r},$i,1) ne "-") {
                $j--;
            }
            $i++;
        }

        #Exclude duplicate indels (potentially from seeded positions that align as a single indel locus)
        foreach my $gip (@{$geneindelpos{$gene}}) {
            if($i >= $gip-$excluderadius && $i <= $gip+$excluderadius) {
                $exclude{$indel}=1;
            }
        }

        if(!$exclude{$indel}) {
            push(@{$geneindelpos{$gene}},$i);	    
            $numgoodindels++;
            
            if($intree) {
                if(scalar(keys %branches) > 0) {
                    $numintree++;
                } else {
                    $numouttree++;
                }
            }

            my $multleft = $i - $viewsize;
            my $multlen = ($viewsize*2) + 15;
            if($multleft < 0) {
                $multleft = 0;
            }
            
            if($minbeta{$indel} <= $beta) {
                for(my $alphacut=$pboundl{$indel}; $alphacut<=$pboundu{$indel}; $alphacut++) {
                    $alphacutcount{$alphacut}++;
                }
            }
            
            my $branchids = "";
            if($intree) {
                $branchids = "".join(",", (sort {$a<=>$b} keys %branches));
            }

            printf ">%s: %d,%d,%d,%d [%s]%s\n", $indel,$pboundl{$indel},$pboundu{$indel},$minbeta{$indel},$minscore{$indel},($numgroups{$indel}==2 ? $indelgapsize{$indel} : "?"), ($intree ? " [".$branchids."]" : "");

            foreach my $r (sort {$indelprotgroup{$indel.":".$a}<=>$indelprotgroup{$indel.":".$b}} keys %refstr) {
                my $groupinfo = "";
                if($indelprotgroup{$indel.":".$r}) {
                    $groupinfo = sprintf "%s", $indelprotgroup{$indel.":".$r};
                }
                printf "%".$maxsplen."s\t%16s\t%s %s\n", $spname{$protsp{$r}}, $r, substr($refstr{$r}, $multleft, $multlen), $groupinfo;
            }
            
            if($intree) {
                if(scalar(keys %branches) < 1) {
                    $branchcountunique{"-1"}++;
                    $branchcount{"-1"}++;
                }
                foreach my $b (keys %branches) {
                    if(scalar(keys %branches) == 1) {
                        $branchcountunique{$b}++;
                    }
                    $branchcount{$b}++;
                }
            }


            #Add good indels to phylip strings
            if($indelgapsize{$indel} >= $minindellen &&
               $pboundl{$indel} <= $alpha &&
               $pboundu{$indel} >= $alpha && 
               $minbeta{$indel} <= $beta &&
               $numtaxa{$indel} >= $mintaxa) {

                $goodgene_ph{$gene} = 1;

                if($intree) {
                    if(scalar(keys %branches) < 1) {
                        $branchcountunique_good{"-1"}++;
                        $branchcount_good{"-1"}++;
                    }
                    foreach my $b (keys %branches) {
                        if(scalar(keys %branches) == 1) {
                            $branchcountunique_good{$b}++;
                        }
                        $branchcount_good{$b}++;
                    }
                }

                if($intree && scalar(keys %branches) > 0) {
                    $numintree_ph++;
                } else {
                    $numouttree_ph++;
                }

                $numphchrs++;
                my $cursplitstr = "";
                my $lastchr = 0;
                my %mapchr = ();
                $mapchr{"?"} = "?";
                
                foreach my $s (sort keys %allspecies) {
                    my $foundspr = 0;
                    my $groupchr = "";
                    foreach my $r (keys %refstr) {
                        if(!$foundspr && $protsp{$r} eq $s && $indelprotgroup{$indel.":".$r} > 0) {
                            $foundspr = 1;
                            $groupchr = $indelprotgroup{$indel.":".$r};
                            if($groupchr > 9) {
                                $groupchr = chr(65 + $groupchr - 10);
                            }
                            if(!((0 + $groupchr) > 0)) {
                                $groupchr = "?";
                            }
                            if(length("".$groupchr)!=1) {
                                die "Group chr not size 1: \"$groupchr\"\n";
                            }
                        }
                    }
                    if(!$foundspr) {
                        $groupchr = "?";
                    }
                    if(!($groupchr eq "?") && !$mapchr{$groupchr}) {
                        $lastchr++;
                        $mapchr{$groupchr} = $lastchr;
                    }
                    $phstr{$s} .= $mapchr{$groupchr};
                    $cursplitstr .= $mapchr{$groupchr};
                }

                $splitstr{$cursplitstr}++;
            }
        }
    }
}

my $excludedindels = scalar(keys %exclude);
my $numgenes = scalar(keys %allgenes);
my $numgoodgenes = scalar(keys %goodgenes);
my $numgoodgenes_ph = scalar(keys %goodgene_ph);

print "\n#Genes: $numgoodgenes/$numgenes, Indels: $numgoodindels/$numindels (Repeats: $excludedindels)";
if($intree) {
    print ", Fit Tree: $numintree/".($numouttree+$numintree);
}
print "\n#Genes: $numgoodgenes_ph/$numgenes, Indels: $numphchrs/$numindels (Gap size >= $minindellen, Alpha = $alpha, Beta = $beta)";
if($intree) {
    print ", Fit Tree: $numintree_ph/".($numouttree_ph+$numintree_ph);
}

print "\n#Phylip Discrete Characters:\n";
print "".scalar(keys %allspecies)." $numphchrs\n";

foreach my $s (sort keys %allspecies) {
    my $spabv = $spname{$s}."__________";
    printf "%.10s %s\n", $spabv, $phstr{$s};
}

foreach my $ss (keys %splitstr) {
    unless($ss =~ /\?/) {
        foreach my $ss2 (keys %splitstr) {
            my $badmatch=0;
            for(my $i=0; $i<length($ss); $i++) {
                if(!(substr($ss2,$i,1) eq "?") && !(substr($ss,$i,1) eq substr($ss2,$i,1))) {
                    $badmatch=1;
                }
            }
            if(!$badmatch) {
                $splitstrincl{$ss}+=$splitstr{$ss2};
            }
        }
    }
}

if($intree) {
    print "\nBranch Count: Branch\tExcl (Incl)\tPhylip Excl (Incl)\n";
    for(my $b=-1; $b<$numgroupings; $b++) {
        printf "%d\t%d (%d)\t%d (%d)\n", $b, 0+$branchcountunique{$b}, 0+$branchcount{$b}, 0+$branchcountunique_good{$b}, 0+$branchcount_good{$b};
    }
}

print "\nSplits:\n";
foreach my $ss (sort {$splitstr{$b}<=>$splitstr{$a}} keys %splitstr) {
    if(!($ss =~ /\?/)) {
        printf "%s %d (%d)\n", $ss, $splitstr{$ss}, $splitstrincl{$ss};
    }
}
print "\n";
foreach my $ss (sort {$splitstr{$b}<=>$splitstr{$a}} keys %splitstr) {
    if($ss =~ /\?/) {
        printf "%s %d\n", $ss, $splitstr{$ss};
    }
}

print "\nAlpha Cutoffs (Beta<=$beta):\n";
foreach my $ac (sort {$a<=>$b} keys %alphacutcount) {
    printf "%d %d\n", $ac, $alphacutcount{$ac};
}

eraseprogress();
