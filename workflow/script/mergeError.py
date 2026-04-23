

import os

import argparse

def merge_overlapping_regions(regions):
    merged_regions = []
    if regions == []:
        return merged_regions
    
    sorted_regions = sorted(regions, key=lambda x: x[0])
    init_region = sorted_regions[0]
    
    for i in range(len(sorted_regions) - 1):
        if sorted_regions[i + 1][0] <= init_region[1]:
            init_region[1] = max(init_region[1], sorted_regions[i + 1][1])
        else:
            merged_regions.append(init_region)
            init_region = sorted_regions[i + 1]
    merged_regions.append(init_region)
    
    return merged_regions

def getContigLength(fai_file):
    contig_length = {}
        
    with open(fai_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            contig_length[items[0]] = int(items[1])

    return contig_length

def main():
    parser = argparse.ArgumentParser(description="Merge Errors from NucFlag")
    parser.add_argument("-hm", "--hifi_misassemblies", help="", required=True)
    parser.add_argument("-om", "--ont_misassemblies", help="", required=True)
    parser.add_argument("-f", "--filter_types", help="", default="")
    parser.add_argument("-b", "--coordinates_file", help="", default="")
    parser.add_argument("-m", "--merge_dis", help="", type=int, default=500000)
    parser.add_argument("-r", "--ref_file", help="", required=True)
    parser.add_argument("-o", "--outdir", help="", required=True)
    parser.add_argument("-of", "--outfile", help="", required=True)
    parser.add_argument("-or", "--out_error_region_file", help="", required=True)
    
    args = parser.parse_args()
    
    hifi_misassemblies = args.hifi_misassemblies
    ont_misassemblies = args.ont_misassemblies
    filter_types = args.filter_types.split(',')
    coordinates_file = args.coordinates_file
    ref_file = args.ref_file
    merge_dis = args.merge_dis
    
    outdir = args.outdir
    outfile = args.outfile
    out_error_region_file = args.out_error_region_file
    
    all_error_regions = {}
    with open(hifi_misassemblies,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            if line.startswith('#'):
                continue
            items = line.split('\t')
            error_type = items[3]
            if error_type in filter_types:
                continue
            
            contig = items[0]
            if contig not in all_error_regions.keys():
                all_error_regions[contig] = []
            
            error_start = int(items[1])
            error_end = int(items[2])
            all_error_regions[contig].append([error_start,error_end])
            
    if ont_misassemblies != "":
        with open(ont_misassemblies,'r') as f:
            while True:
                line = f.readline()[:-1]
                if not line:
                    break
                if line.startswith('#'):
                    continue
                items = line.split('\t')
                error_type = items[3]
                if error_type in filter_types:
                    continue
                contig = items[0]
                if contig not in all_error_regions.keys():
                    all_error_regions[contig] = []
                error_start = int(items[1])
                error_end = int(items[2])
                all_error_regions[contig].append([error_start,error_end])
    
    outfile = open(outfile,'w')
    
    all_merged_regions = {}
    for i in all_error_regions.keys():
        all_merged_regions[i] = merge_overlapping_regions(all_error_regions[i])

    for i in all_merged_regions.keys():
        for j in all_merged_regions[i]:
            outfile.write(i + '\t' + str(j[0]) + '\t' + str(j[1]) + '\n')
    outfile.close()
    
    # 考虑两种情况，如果coordinates_file提供，如果coordinates_file不提供，如果不提供，则直接将所有错误合并输出，如果提供，则只输出坐标文件中指定区间内的错误
    contigs_regions = {}
    if coordinates_file == "":
        cmd = 'samtools faidx ' + ref_file
        os.system(cmd)
        fai_file = ref_file + '.fai'
        contig_length = getContigLength(fai_file)
                    
        chr_merge_bed = {}
        for i in all_error_regions.keys():
            chr_merge_bed[i] = []
            if len(all_error_regions[i]) == 1:
                chr_merge_bed[i]  = all_error_regions[i]
            else:
                sorted_bed = sorted(all_error_regions[i],key=lambda x:x[0])
                # 先去掉重复,如果下一个的end小于上一个，则删掉
                remove_dup = []
                remove_dup.append(sorted_bed[0])
                for j in range(len(sorted_bed) - 1):
                    if sorted_bed[j + 1][1] <= sorted_bed[j][1]:
                        continue
                    else:
                        remove_dup.append(sorted_bed[j + 1])
                # 合并
                init_region = remove_dup[0]
                for j in range(len(remove_dup) - 1):
                    if (remove_dup[j + 1][0] - init_region[1]) < merge_dis:
                        init_region = [init_region[0], remove_dup[j+1][1]]
                    else:
                        chr_merge_bed[i].append(init_region)
                        init_region = remove_dup[j + 1]
                chr_merge_bed[i].append(init_region)
        final_regions = []
        for i in chr_merge_bed.keys():
            for j in chr_merge_bed[i]:
                start = int(j[0]) - int(merge_dis / 2)
                if start < 1:
                    start = 1
                end = int(j[1]) + int(merge_dis / 2)
                if end > contig_length[i]:
                    end = contig_length[i]
                final_regions.append([i, start, end])
                if i not in contigs_regions.keys():
                    contigs_regions[i] = {}
                if (start,end) not in contigs_regions[i].keys():
                    contigs_regions[i][(start,end) ] = []
                contigs_regions[i][(start,end)] = []        
    else:
        # 如果提供了坐标文件，则只处理坐标文件中指定区间内的错误
        with open(coordinates_file,'r') as f:
            while True:
                line = f.readline()[:-1]
                if not line:
                    break
                items = line.split('\t')
                contig = items[0]
                region_start = int(items[1]) - 1
                region_end = int(items[2])
                if contig not in contigs_regions.keys():
                    contigs_regions[contig] = {}
                if (region_start,region_end) not in contigs_regions[contig].keys():
                    contigs_regions[contig][(region_start,region_end)] = []
    
    for i in all_merged_regions.keys():
        for j in all_merged_regions[i]:
            contig = i
            error_start = j[0]
            error_end = j[1]
            for region in contigs_regions[contig].keys():
                region_start = region[0]
                region_end = region[1]
                if error_start >= region_start and error_end <= region_end:
                    contigs_regions[contig][region].append([error_start - region_start, error_end - region_start])
                    break
    out_error_region_file = open(out_error_region_file,'w') 
    for i in contigs_regions.keys():
        for j in contigs_regions[i].keys():
            if len(contigs_regions[i][j]) == 0:
                continue   
            out_error_region_file.write(i + '\t' + str(j[0] + 1) + '\t' + str(j[1]) + '\n') 
            outfile = outdir + '/' + i + ':' + str(j[0] + 1) + '-' + str(j[1]) + '.xls'
            outfile = open(outfile,'w')
            for k in contigs_regions[i][j]:
                outfile.write(str(k[0]) + '\t' + str(k[1]) + '\n')
            outfile.close()
            
    out_error_region_file.close()
           
    

if __name__ == "__main__":
    main()