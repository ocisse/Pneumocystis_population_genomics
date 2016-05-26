# need to comment this better

# produce clusters and output and align

ls *.merged_plus_outgroup.reordered.aln > list.txt

# produce stats
perl Polymorphorama_coding.pl list.txt frquency_cutoff 1
