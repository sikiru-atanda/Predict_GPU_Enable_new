plink2 --vcf FFAR_Chinese_Pea_ref_genome_filtered.vcf.gz --allow-extra-chr \
        --min-alleles 2 --max-alleles 2 --geno 0.10 --maf 0.05 --export vcf bgz --out initial_qc_data
