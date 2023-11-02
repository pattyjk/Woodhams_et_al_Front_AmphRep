## Cluster against 2023 version of Amphibac Database
```
#export rep seqs from QIIME2
conda activate qiime2-2023.5
qiime tools export --input-path Filt_Aligned_Repset_seqs.qza --output-path rep_seqs
conda deactivate

#cluster at 100% against database
conda activate vsearch
vsearch -usearch_global rep_seqs/rep_seqs.fna -db antiBd_db/AmphibBac_Inhibitory_2023.2r.fasta --strand plus --id 0.97 --blast6out leopard_frog_out.txt

#get same results (18.42% of sOTUs) match the database between 0.91 and 1 percent identity. 
```