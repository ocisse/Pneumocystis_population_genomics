#!/nethome/cisseoh/perl5/perlbrew/perls/perl-5.20.3-thread-multi/bin/perl -w
#
#
use strict;
use Data::Dumper;
use Carp;
use IO::All;
use feature 'say';

my @files = load($ARGV[0]);

my %moiData = ();
my $i = ""; 
foreach $i (@files){
	my @data = extract_moi($i);
	@{$moiData{$i}} = @data
}

say "FILE,MOI,Count,Total";
my $o = "";
foreach $o (keys %moiData){
	my @tmp = @{$moiData{$o}};
	print "$o";
	my $x = "";
	
	if ($tmp[0]){
		foreach $x (@tmp){
			print ",$x";
		}
		print "\n";
	} else {
			print ",0,0,0\n";
	}
}
# sub
sub extract_moi {
	my @data1 = ();
	my ($file) = @_;
	my $f1 = io($file);
           $f1->autoclose(0);
           while(my $l1 = $f1->getline || $f1->getline){
            chomp $l1;
		if ($l1 =~ m/MOI-estimate/){
			my @d = split /\t/, $l1;
			my ($moi,$count,$total) = @d[0..2];
			push(@data1,$moi,$count,$total);		
		} else {
			next;
		}
        }
	return(@data1)
}
sub load {
	my @ar = ();
	my $f = io(@_);
	   $f->autoclose(0);
	    while(my $l = $f->getline || $f->getline){
	    chomp $l;
		push(@ar,$l);
	}
	return(@ar);
}
