# Whole genome analysis of Histiotus bats
These are the programs, commands and files used in the paper "Whole-genome analysis reveals contrasting relationships among nuclear and mitochondrial genomes between three sympatric bat species" by Laine et al 2023 https://doi.org/10.1093/gbe/evac175

Sequence data can be found in Genbank under Bioproject PRJNA827658.  
Assembled mitochondria are in Genbank under accession number OP328298-OP328300 and OP345288-OP345340.

1. Read mapping 

Program used bwa v. 0.7.1 against Eptesicus fuscus genome (GCA_000308155.1 EptFus1.0): bwa mem ref.fa reads.fq -B 3 -O 5 -k 15

2. Genotype likelihood calling 

Program used ANGSD 0.935: angsd -GL 1 -out GCA_000308155.1_EptFus1.0_genomic.fna -doGlf 2 -doMajorMinor 1 -doMaf 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -setMaxDepth 100 -minMapQ 20 -minQ 20 -minMaf 0.05 -minInd 30[only histiotus]/40[all] -doCounts 1 -bam filelist_histi.txt

3. LD pruning 

Program used ngsLD 1.1. with prune_graph.pl using parameters max_kb_dist 5 --min_weight 0.5.

4. PCA 

Program used PCAngsd: pcangsd.py -beagle pruned_final.beagle.gz -out PCA

5. Admixture 

Program used NgsAdmix v.32:  
for k in 'seq 1 3'; do for eachK in 'seq 1 10';  
do  
NGSadmix -likes pruned_final.beagle.gz -K "$k" -P 8 -o admix_"$k"_"$eachK" -minMaf 0.05  
done; done  
Following CLUMPAK http://clumpak.tau.ac.il/bestK.html for best K.  
Good tutorial: https://github.com/alexkrohn/AmargosaVoleTutorials/blob/master/ngsAdmix_tutorial.md

6. Genetic diversity 

- Make SAF -file per Histiotus species with ANGSD:  
angsd -b hmac_filelist.bams -ref GCA_000308155.1_EptFus1.0_genomic.fna -anc GCA_000308155.1_EptFus1.0_genomic.fna -out Hmac_saf -P 10 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 10 -setMaxDepth 100 -doCounts 1 -GL 1 -doSaf 1

- 2dSFS per Histiotus species pair with ANGSD:  
realSFS -P 20 Hmac_saf.saf.idx Hmag_saf.saf.idx -nSites 500000000 > Hmac_Hmag_500.sfs

- Pairwise Fst with ANGSD:  
realSFS fst index Hmac_saf.saf.idx Hmag_saf.saf.idx -sfs rowsum_Hmac_Hmag_500.sfs -fstout Hmac_Hmag_500.fst -whichFst 1  
realSFS fst stats Hmac_Hmag_500.fst.fst.idx

- 1dSFS and diversity indices per species with ANGSD:  
realSFS -P 20 Hmac_saf.saf.idx -nSites 500000000 > Hmac_500.sfs  
realSFS saf2theta Hmac_saf.saf.idx  -outname Hmac -sfs Rowsum_Hmac_500.sfs  
thetaStat do_stat Hmac.thetas.idx -outnames Hmac_stats  

7. Genetic introgression 

Patterson’s D statistics were calculated in ANGSD using the ABBABABBA2 method: angsd -b bamlist -ref GCA_000308155.1_EptFus1.0_genomic.fna -doMajorMinor 1 -GL 1 -doMaf 1 -doCounts 1 -doAbbababa2 1 -sizeFile pop.size -useLast 1 -rf regions_1000scaffs.txt -baq 1 -remove_bads 1 -uniqueOnly 1 -C 50 -minMapQ 20 -only_proper_pairs 1 -minQ 20 -minInd 30 -setMaxDepth 100 -SNP_pval 1e-6 -P 10 -out AB_histiebot_pop

The D-values were called with the R-script estAvgError provided by ANGSD.  

8. Nuclear phylogenetic tree 

A neighbor-joining tree was constructed with the BioNJ tree building algorithm of FastME v.2.1.5 (Lefort et al. 2015), based on individual pairwise genetic distances estimated with ngsDist v.1.0.9: ngsDist --posH all_beagle_pos.txt --geno genolike_all.beagle.gz --probs --n_ind 58 --n_sites 9943522 --n_boot_rep 100 --n_threads 10 --out all_dist_boot100  
FastME http://www.atgc-montpellier.fr/fastme/ with BioNJ, input provided from previous step.  
Support values with RAxML-NG v. 1.0.2: raxml-ng --support --tree all_dist_boot100.main.nwk --bs-trees all_dist.boots100.nwk --prefix all_dist_boot100  

A maximum likelihood tree was constructed with a SNP panel. SNPs were called from the genotype likelihood file created above (all species) with BEAGLE Utilities https://faculty.washington.edu/browning/beagle_utilities/utilities.html program gprobs2beagle using minimum posterior probability of 0.8. Then beagle2vcf from BEAGLE Utilities was used to transform the file to a vcf-format which was further transformed to PHYLIP-format with vcf2phylip v. 2.0 (Ortiz 2019) with minimum sample locus of 20. This provided 8 511 209 SNPs. Tree was constructed with IQ-TREE v. 2.1.4_beta (Minh et al. 2020) with model finder (Kalyaanamoorthy et al. 2017) and bootstrapping (1000) (-bb 1000 -m TEST). 

9. Mitochondrial DNA assembly per species 

Whole mitochondria were assembled with GetOrganelle v. 1.7.5.: get_organelle_from_reads.py -1 R1.fastq.gz -2 R2.fastq.gz -R 10 -k 21,45,65,85,105 -F animal_mt -o Sample_1  

If some samples fail to assemble, increase last k to 127 (-k 21,45,65,85,127) and/or add seed (-s good_sample.fasta).

10. Mitochondrial phylogeny 

Assembled mitochondria and E. fuscus 215 mitochondria from NCBI (MF143474.1) were aligned with Clustal Omega v. 1.2.4: clustalo -i all_mt_rehead.fasta -o all_mt_clustalo.fasta --distmat-out=all_mt_clustalo.dist --full

The Clustal alignment file was separated into 13 protein-coding mitochondrial gene and 2 rRNA alignments based on the E. fuscus mitochondrial annotation and combined in one Nexus file with help of Geneious 11.03. Overlaps were also removed between ATP8 and ATP6, ND4L and ND4, ND5 and ND6 to include the overlapping region only once.

A consensus tree was built with IQ-TREE v. 2.1.4_beta, input provided: iqtree -T AUTO -s All_bygene_clustalo_TRIM_rRNAs.nex -p charset_TRIM_rRNAs.nex -bb 1000 -m MFP+MERGE -pre all_mtgenes_part_TRIM_rRNAs_iqtree

11. Mitochondrial molecular timing 

BEAST2 v. 2.6.6. for molecular dating and BEAUti 2 to produce the run file for BEAST2. BEAST2 was ran three times with chain length 50 000 000 and the trees were combined with LogCombiner with 10% burnins and TreeAnnotator was used for consensus tree. Input file (all_bygenes_partition_50m_inputREVISION.xml) is provided. 
