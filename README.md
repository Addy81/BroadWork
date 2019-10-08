
# What have I done?

VCF comparison script:
* Simultaneously reads and compares two vcf files from the same sample
* **Required packages:**
  * PyVCF (pip install PyVCF)
  * numpy (pip install numpy)
* **Input**: two vcf files to be compared (should work for zipped and unzipped)
  * If the sample name is different than the vcf sample name, provide it as a third argument.
* **Returns**: 
  * List of vcf entries that have different DP values
  * Summary of how many records in the vcf have different DP values and different records all together.
  * Histogram of DP difference variance
  * List of different vcf records (record = [CHROM,POS,REF,ALT])
*  **Runtime**: can take up to 5 hours for genome vcfs depending on the size of the vcf.
   *  If so the vcfs can be split by chromosome


# How to - step by step

## Run the script

    $ python vcf_compare.py original_vcf.gz development_vcf.py > output_filename.txt

The tool gets the sample name from the filename, it has to be the same sample name as it is in the vcf, if it isn't, provide a sample name as a third argument.

    $ python vcf_compare.py original_vcf.gz development_vcf.py NA12878 > output_filename.txt

### For big files:

   1. Index your vcf files using tabix ()
   
    $ tabix -p vcf your_vcf_file.vcf.gz
   2. Create a file listing all chromosomes set it as a bash variable.

    $ chromosomes=/path/to/file/chromosomes.txt
   3. Split the vcfs by chromosome
   
    $ for chrom in $(cat $chromosomes); do echo ${chrom} splitting; tabix -h NA12878.truth.g.vcf.gz ${chrom} > NA12878.dev.${chrom}.g.vcf; done

   4. Optional : recompress

    $ gzip *.vcf

5.  Run script on indicidual components

```
$ for chrom in $(cat $chromosomes); do echo $chrom running; python ${path}/bin/vcf_compare.py NA12878.truth.${chrom}.g.vcf.gz NA12878.dev.${chrom}.g.vcf.gz > NA12878.${chrom}.comparison.txt ; done
```
