#!/nethome/cisseoh/perl5/perlbrew/perls/perl-5.20.3-thread-multi/bin/perl -w
#
use strict;
use Data::Dumper;
use IO::All;
use feature 'say'; 
use Carp;


my $f = io($ARGV[0]); 
   $f->autoclose(0);
   while(my $l = $f->getline || $f->getline){
   chomp $l;
	next if $l =~ m/^#/;
	my @d = split /\t/,$l;
	
	my ($wStart,$wEnd) = $d[0] =~/\)\((\d+)\,(\d+)\)\(/;
	say "$d[1]\t$wStart\t$wEnd\t$d[8]";

}
