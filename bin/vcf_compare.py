#!/usr/bin/python3
#
#
#
#
# Adriana Toutoudaki (September 2019) contact: atoutoud@broadinstitute.org

import vcf
import os, sys
import numpy as np
from vcf.utils import walk_together

class Comparison:

    def __init__(self):
        self.matches = 0
        self.diff_metrics = 0
        self.diff = 0
        self.no_format_count = 0
        self.no_formats = []
    
    def total_count(self):
        """ Returns total count of calls per vcf."""
        self.tc = self.matches + self.diff + self.no_format_count +self.diff_metrics

        return self.tc

    def get_stats(self):
        """
        Calculates the percentage of calls with different depth and different records.
        Returns those values.
        """
        depth_stats = (self.diff_metrics / self.tc) * 100
        diff_stats = (self.diff/ self.tc) * 100

        return depth_stats,diff_stats

    def output_no_format(self):
        """ Currently not used in the code below"""
        
        with open('reads_missing_format_field.txt','w+') as f:
            header_line = '\t'.join(['Chrom','Pos','Ref','Alt','INFO'])
            f.write(header_line)

            for read in self.no_formats:
                line = '\t'.join([read[0].CHROM,str(read[0].POS),read[0].REF,read[0].ALT])
                f.write(line)

    def print_metrics(self,depth_stats,diff_stats):
        """ Prints summary of differences."""
        print ('Total\tDiffDP\tDiffRec')
        print('{}\t{}\t{}'.format(self.tc,self.diff_metrics,self.diff))
        print('{}\t{}\t{}'.format('100%',depth_stats,diff_stats))

def calculate_hist(lst,bin_edges):
    """
    param1: list of DP differences or DP % change gathered during comparison
    param2: desired bin edges, for DP difference histogram [1,2,3,4,5,6,7,8,9,max(DP_difference_list)] 
    for % difference histogram [0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,max(percent_difference_list)])
    prints: Histogram, bin edges and % off calls in each bin. 
    """
    abs_lst = [abs(x) for x in lst]
    counts, bin_edges = np.histogram(abs_lst,bins=bin_edges)
    per_lst = [round(count*100/sum(counts),2) for count in counts]
    
    #Print the outliers to emphasise the contents of the last bin
    outliers = []
    for x in abs_lst:
        if x >8:
            outliers.append(x)
    
    print ('\nHistogram:', counts)
    print ('Bins:',bin_edges)
    print ('Percentage of values in each bin:', per_lst)
    
    print ('Values in the last bin:',outliers)

def main():


    #dev_vcf = '/Users/atoutoud/Projects/testCompare/data/NA12891.dev_short.vcf'
    #truth_vcf = '/Users/atoutoud/Projects/testCompare/data/NA12891.truth_short.vcf'

    truth_vcf = sys.argv[1]
    dev_vcf = sys.argv[2]

    #dev_reader= vcf.Reader(open(dev_vcf, 'r'))
    #truth_reader = vcf.Reader(open(truth_vcf,'r'))

    print ('Sample Comparison')
    dev_reader= vcf.Reader(filename=dev_vcf)
    truth_reader = vcf.Reader(filename=truth_vcf)

    #Checks if a filename is provided. pyVCF looks for the filename in the header line, for replicates of the same sample with different filenames
    # the correct ones should be provided otherwise it will fail.
    if len(sys.argv) == 4:
        sample = sys.argv[3]
    else:
        sample = os.path.basename(truth_vcf).split(os.extsep)[0]
    
    print (sample)
    summary = Comparison()
    records_dont_match = []
    call_difference = []
    percent_difference = []
    DP_range = []

    #Walk_together is a pyVCF inbuilt function to read two vcfs at the same time.
    for dev_rec,truth_rec in walk_together(dev_reader,truth_reader):
        # A record corresponds to [CHROM,POS,REF,ALT], if the same it checks the metrics differences.
        if dev_rec == truth_rec:
            try:
                #If the DP is different between the records.
                if dev_rec.genotype(sample)['DP'] != truth_rec.genotype(sample)['DP']:
                    summary.diff_metrics +=1 
                    #count_metrics += 1
                    print('')
                    print (dev_rec.CHROM,dev_rec.POS, dev_rec.REF,dev_rec.ALT,dev_rec.QUAL)
                    print ('--------------------------------------------------------------')
                    print ('\t'.join(dev_rec.FORMAT.split(':')))
                    
                    for entry in truth_rec.genotype(sample).data:
                        print (entry,end = '\t')
                    print('')
                    
                    for entry in dev_rec.genotype(sample).data:
                        print (entry,end='\t')
                    
                    true_DP = truth_rec.genotype(sample)['DP']
                    test_DP = dev_rec.genotype(sample)['DP']
                    DP_range.append(true_DP)
                    DP_range.append(test_DP)

                    if true_DP == 0:
                        difference = 0 # had to set, as the % different calculation divides by the true_DP, so if that is 0 it breaks
                    else:
                        difference = round(abs((test_DP-true_DP)/true_DP*100),4)
                    
                    percent_difference.append(difference)

                    DP_diff = true_DP-test_DP
                    call_difference.append(DP_diff)
                    
                    print ('\nDP difference {}'.format(DP_diff))
                    print (difference,'%')
                    print ('')

                else: 
                    summary.matches +=1
            
            #
            except AttributeError: 
                summary.no_format_count +=1
                summary.no_formats.append([dev_rec,dev_rec.INFO])
                print ('No format fields {} at position:{}'.format(dev_rec.CHROM,dev_rec.POS))
            
        else:
            #count_no_match +=1
            summary.diff +=1
            
            #Stores the different values so they can be explorted all together at the end. 
            if truth_rec is None:
                records_dont_match.append({"truth":(truth_rec),"dev":(dev_rec.CHROM,dev_rec.POS,dev_rec.REF,dev_rec.ALT)})
            elif dev_rec is None:
                records_dont_match.append({"truth":(truth_rec.CHROM,truth_rec.POS,truth_rec.REF,truth_rec.ALT),"dev":(dev_rec)})
            else:
                records_dont_match.append({"truth":(truth_rec.CHROM,truth_rec.POS,truth_rec.REF,truth_rec.ALT),"dev":(dev_rec.CHROM,dev_rec.POS,dev_rec.REF,dev_rec.ALT)})
                 #print ('** Records do not match **',dev_rec,truth_rec)


    summary.total_count()

    stats1,stats2 = summary.get_stats()
    #summary.output_no_format()
    
    #Prints the summary metrics for the entirety of the vcf files
    summary.print_metrics(round(stats1,4),round(stats2,4))

    print ("\nRecords that didn't match first is truth, second is dev")
    for i in records_dont_match:
        print (i)

    #Outputs histogram values for DP difference and percent difference.
    calculate_hist(call_difference,[1,2,3,4,5,6,7,8,9,max(call_difference)])
    print ('***Please note last bin contains values of difference greater than 8 calls***')
    calculate_hist(percent_difference,[0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,max(percent_difference)])
    print ('***Please note last bin contains entities with a percent change greater than 5% ***')


    print ('\nDP range:',min(DP_range),max(DP_range))
if __name__ == '__main__':
    main()
