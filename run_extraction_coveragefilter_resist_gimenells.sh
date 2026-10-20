

echo "Running freebayes"

freebayes -f dmel-all-chromosome-r6.12.fasta \
 --pooled-continuous \
 --targets resistance_genes.bed \
 --min-coverage 246 \
 -g 1477 \
 -L list.txt \
 > variants_allc_pooled_cov_resist.vcf


# Komprimiert die grosse VCF
echo "Compressing VCF-file"
bgzip -f variants_allc_pooled_cov_resist.vcf


# Index
echo "Creating index for compressed VCF-file"
tabix -f variants_allc_pooled_cov_resist.vcf.gz


# Chromosom, Position, Ref. Allel, altern. Allele, allele depths (Reads)
echo "Extracting allele counts with bcftools"

bcftools query \
 -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' \
 variants_allc_pooled_cov_resist.vcf.gz \
 > allele_counts_snps_cov_resist.tsv


echo "Filtering with awk"

awk 'BEGIN { FS="\t"; OFS="\t" }
{
    # Split ALT field (column 4) on commas.
    n = split($4, altAlleles, ",");

    # Only consider sites with at least 2 alternate alleles
    # (REF + >=2 ALTs = at least 3 alleles)
    if(n >= 2) {
        valid = 1;

        # Check that all sample fields are not missing
        for(i = 5; i <= NF; i++){
            if($i == ".") {
                valid = 0;
                break;
            }
        }

        if(valid)
            print $0;
    }
}' allele_counts_snps_cov_resist.tsv \
> filtered_sites_cov_resistence_genes.tsv


echo "Done!"