#!/sysapps/cluster/software/Perl/5.20.2-goolf-1.7.20/bin/perl  -w

use strict; 
use Data::Dumper; 
use feature 'say'; 
use IO::All; 
use Carp;

# my refernece CDS, 

my $pjRefCDS = 'GCA_001477535.1_Pneu_jiro_RU7_V2_genomic.CDS.fa'; 
my $pcRefCDS = 'GCA_001477545.1_Pneu_cari_B80_V3_genomic.CDS.fa';
my $pmRefCDS = 'GCF_000349005.1_Pneumo_murina_B123_V2_genomic.CDS.fa';

# output
my $pj_Pc_outgroup_dir = 'Pj_Pc_outgroup_dir';

# load orthology
my %og = load_orthologs("out");
#say Dumper \%og;

# load files I need to process
my @pjfiles = loadFiles('Pj_files_list.txt'); # to update
#say @pcfiles;

# first buil Pc - Pm outgroup
my @pjRefIds = extractRefIds($pjRefCDS);
my $i = "";
foreach $i ( @pjRefIds ) {
	my ($cleanId) = $i =~/(T551_\d+)/;
	my %problematic_ids = load_excluded_ids('Pj_excluded_locus_tags.txt');
	
	next if (defined($problematic_ids{$cleanId}));

	my $orhtolog = `grep $cleanId out`;
	if ( $orhtolog ne '') {
		$orhtolog =~s/\n//;
		my @og = split /\t/, $orhtolog;
		my ($pc,$pj,$pm) = @og[0..2];

		$pc =~s/T0//;
		$pj =~s/T0//;
		$pm =~s/T0//;
		
		if ( -e "$pj_Pc_outgroup_dir/$cleanId.merged_plus_outgroup.fa") {
			warn"$pj_Pc_outgroup_dir/$cleanId.merged_plus_outgroup.fa\t exists\t...skip\n\n";
		} else {
			write_this_locus_tag_data($cleanId,'Pj_files_list.txt',$pjRefCDS, $pc, $pcRefCDS, $pj_Pc_outgroup_dir);
		}
	}
}

# sub
sub load_excluded_ids{
	my %he = ();
	my $e = io(@_);
	   $e->autoclose(0);
	   while (my $el = $e->getline || $e->getline ){
	   chomp $el;
		$he{$el} = $el;
	}
	return(%he);
}  
sub write_this_locus_tag_data {
	my ($locus_tag,$listOfFiles,$ref,$outgroupid,$outgroupReffile, $outputdir) = @_; 

	# extract the locus from each fasta # key = abbre for sample names vla = name of the files
	my %files = loadFiles($listOfFiles); 

	my $s = "";
	foreach $s ( keys %files ) {
	#	my $id = `grep $locus_tag $s`;
		my $id = extractHeader($locus_tag,$s);	
		
		if (defined($id)) {
			run("echo $id > $outputdir/$s-$locus_tag.tmp");
			run("perl ~/utils/retrieveFasta.pl $outputdir/$s-$locus_tag.tmp $s > $outputdir/$s-$locus_tag.fa");
		
			# trying to place the name of the sample inside 
			my @cleanid = split /_/,$id;
			my $cleanIdtmp =  "$cleanid[0]_$cleanid[1]";
			my $cleanheader	 = "PJIROVECII_$files{$s}_$cleanIdtmp";  
			run("perl -pi.old -E \"s/$id/$cleanheader/\" $outputdir/$s-$locus_tag.fa") if (defined($id));;
		}
	}	
	# extract the outgroup sequence
	my $outg = extractHeader($outgroupid,$outgroupReffile);
	my @outtmp = split /\s/, $outg;
	run("echo $outtmp[0] > $outputdir/$outtmp[0].tmp");
	run("perl ~/utils/retrieveFasta.pl $outputdir/$outtmp[0].tmp $outgroupReffile > $outputdir/$outtmp[0].fa");	
	
	# concat    
	run("cat $outputdir/$outtmp[0].fa $outputdir/*-$locus_tag.fa > $outputdir/$locus_tag.merged_plus_outgroup.fa");


	# align
	run("muscle -in $outputdir/$locus_tag.merged_plus_outgroup.fa -out $outputdir/$locus_tag.merged_plus_outgroup.aln");	
	
	# need to re-order - the outsequence need to be the first sequence in the file
	# remove sequence that contains 50% Ns
	check_AND_reorder_fasta("$outputdir/$locus_tag.merged_plus_outgroup.aln","$outputdir/$locus_tag.merged_plus_outgroup.reordered.aln");

	# clean
	run("rm $outputdir/*-$locus_tag.* $outputdir/$outtmp[0].* ");
}

sub check_AND_reorder_fasta{
	my ($infas,$reorderfas) = @_;
	
	my $buffer = "";

	my %h2 = ();
	open(F,$infas) || die "couldn't read the file: infas:$!\n";
	local $/ = "\n>"; 
	
	while ( my $seq = <F> ) {
	chomp $seq;
		my ($header) = $seq =~/(>*.+)\n/;
		$header =~s/>//;
		$seq =~s/>*.+\n//;
		$seq =~s/\n//g;
		$h2{$header} = $seq;
	}
	close F;
	
	my $u = "";
	foreach $u ( keys %h2 ) {
		if ( $u =~ m/T552/){
			$buffer .=">$u\n";
			$buffer .="$h2{$u}\n";
		} else {
			next;
		}
	}
	foreach $u ( keys %h2 ) {
                unless ( $u =~ m/T552/){
			
			my @alnbases = split //, $h2{$u};
			my $n = 0;
			my $ns = "";
			foreach $ns ( @alnbases ) {
				$n++ if ( $ns eq 'N' );
			}
			
			unless ($n > (scalar @alnbases / 2)){
                        	$buffer .=">$u\n";
                        	$buffer .="$h2{$u}\n";
                	} 
		}
        }
	io("$reorderfas")->write($buffer);
}
sub extractHeader {
	my ($tofind,$fasta) = @_;
	
	my $fas = io($fasta);
           $fas->autoclose(0);
	   while ( my $fasl = $fas->getline || $fas->getline ) {
	   chomp $fasl;
		if ($fasl ~~ /$tofind/) {
			$fasl =~s/>//;
			return($fasl);
			last;
		} else {
			next;
		}	
	}
}
sub run {
	my ($cmd) = @_;
	say "###\t$cmd:$!\n";
	system($cmd)==0 || carp"cannot run $cmd:$!\n";
}
sub extractRefIds {
	my @ar1 = ();
	my $f1 = io(@_);
	   $f1->autoclose(0);
	   while (my $l1 = $f1->getline || $f1->getline ){
	   chomp $l1;
		if ( $l1 =~ m/^>/) {
			$l1 =~s/>//;
			push(@ar1, $l1);
		} else {
			next;
		}

	}
	return(@ar1);
}
sub loadFiles {
	my %h = ();
	my $f = io(@_);
	   $f->autoclose(0);
	   while ( my $l = $f->getline || $f->getline ){
	   chomp $l;
		my @dat = split /\t/, $l;
		$h{$dat[1]} = $dat[0];
	}
	return(%h);
}
sub load_orthologs {
        my %h = ();
        my $c = 0;

        my $f1 = io(@_);
        $f1->autoclose(0);
        while (my $l1 = $f1->getline || $f1->getline ){
        chomp $l1;
                next if $l1 =~ m/^#/;
                next if $l1 =~ m/\*/; # skip missing in one of the pneumo
                my @tmp = split /\t/, $l1;
		
		my @newTmp = (); 
		my $r = "";
		foreach $r ( @tmp ) {
			$r =~s/T0//;
			push(@newTmp,$r);
		}
                @{$h{"OG.$c"}}  = @newTmp;
                $c++;
        }
        return(%h);
}
