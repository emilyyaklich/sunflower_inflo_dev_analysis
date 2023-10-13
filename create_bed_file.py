# Name: create bedfile
# Author: EY
# Date: May 01 2023 (edited October 13 2023)
# Version: Python 3.9
# Description: create bedfile for bedtools input 'getfasta' for HAshoot gtf file

import pandas as pd
import re
# create bedfile for bedtools input 'getfasta'

# chromosome start end name
chrom= []
start=[]
end=[]
name=[]
with open('/home/emilyyaklich/Documents/sunflower_inflo_dev_analysis/Hashoot_ref_merged.gtf') as topo_file:
    for line in topo_file:
        if '#' not in line:

            split_line = re.split(r'\t+', line)
            if "transcript" in split_line[2]:
                start.append(int(split_line[3]))
                end.append(int(split_line[4]))
                split_quote = split_line[8].split('"')
                mrna_id = split_quote[1]
                name.append(mrna_id)
                chrom_name=split_line[0]
                chrom.append(chrom_name)

d = {'Chrom': chrom, 'Start': start, 'End':end,'Name':name}
df = pd.DataFrame(d)
df.to_csv('getfasta_bed_input.bed', sep='\t',header=False,index=False)
