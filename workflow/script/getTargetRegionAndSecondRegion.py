
import os
import argparse

def run_getTargetArray(ref_file,cen_list_file,outdir,outdir_cen_array):
    ref_asm = ref_file 
    with open(cen_list_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            array_name = items[0] + ':' + items[1] + '-' + items[2]
            cmd = 'samtools faidx ' + ref_asm + ' ' + array_name + ' > ' +  outdir_cen_array + '/' + array_name + '.tmp.fa'
            os.system(cmd)
            cmd = 'seqkit seq -w 0 ' + outdir_cen_array + '/' + array_name + '.tmp.fa' + ' > ' + outdir_cen_array + '/' + array_name + '.fa'
            os.system(cmd)
            cmd = 'rm ' +  outdir_cen_array + '/' + array_name + '.tmp.fa'
            os.system(cmd)
    cmd = 'cat ' + outdir_cen_array + '/*.fa > ' + outdir + '/all.cenarray.fa'
    os.system(cmd)



def getFilterRecord(paf_file,cen_list_file,outfile,MAPQ_thr,break_contig_merge_thr,min_array_thr,overlap_merge_thr): 
    target_arrays = {}
    
    # 构建最初的target_arrays
    with open(cen_list_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            array_name = items[0] + ':' + items[1] + '-' + items[2]
            target_arrays[array_name] = {}
    
    # 从paf 文件中读取，保存到target_arrays
    # key1: array, key2:hifiasm-contig_strand, value: 对应记录
    with open(paf_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            # print(items)
            if int(items[11]) <= MAPQ_thr:
                continue 
            if items[0] in target_arrays.keys():
                if items[5] + '_' +  items[4] not in target_arrays[items[0]].keys():
                    target_arrays[items[0]][items[5] + '_' +  items[4] ] = []
                    target_arrays[items[0]][items[5] + '_' +  items[4] ].append([items[0],items[1],int(items[2]),int(items[3]),'+',items[5],int(items[7]),int(items[8]), items[4],items[11]])
                else:
                    target_arrays[items[0]][items[5] + '_' +  items[4] ].append([items[0],items[1],int(items[2]),int(items[3]),'+',items[5],int(items[7]),int(items[8]), items[4],items[11]])
    
    
    
    # 对于记录按照hifiasm start进行排序
    sorted_target_arrays = {}
    for i in target_arrays.keys():
        sorted_target_arrays[i] = {}
        for j in target_arrays[i].keys():
            sorted_target_arrays[i][j] = sorted(target_arrays[i][j],key=lambda x:x[6])

    # 对于每一个sorted_target_arrays[key1][key2]，合并其中记录基于，大小break_contig_merge_thr
    target_arrays_tmp = {}
    for i in sorted_target_arrays.keys():
        target_arrays_tmp[i] = {}
        for j in sorted_target_arrays[i].keys():
            target_arrays_tmp[i][j] = []
            if len(sorted_target_arrays[i][j]) == 1:
                target_arrays_tmp[i][j] = sorted_target_arrays[i][j]
            else: 
                merge_records = []
                init_record = sorted_target_arrays[i][j][0]
                for k in range(len(sorted_target_arrays[i][j]) -1):
                    next_record = sorted_target_arrays[i][j][k + 1]
                    if abs(next_record[6] - init_record[7]) < break_contig_merge_thr:
                        if next_record[7] > init_record[7]:
                            # 进行合并
                            init_record[7] = next_record[7]
                            if init_record[2] > next_record[2]:
                                init_record[2] = next_record[2]
                            if init_record[3] < next_record[3]:
                                init_record[3] = next_record[3]
                        else:
                            continue
                    else:
                        merge_records.append(init_record)
                        init_record = next_record
                merge_records.append(init_record)

                target_arrays_tmp[i][j] = merge_records
  

    target_arrays = {}
    # 重新生成记录，并去掉小于100k的
    # 重新存成key: array, value: 记录
    for i in target_arrays_tmp.keys():
        target_arrays[i] = []
        for j in target_arrays_tmp[i].keys():
            for k in target_arrays_tmp[i][j]:
                if k[7] - k[6] > min_array_thr:
                    target_arrays[i].append(k)

    # second contig必须唯一对应target array,如果不唯一，尝试选择比对长度最长的 
    # note: *着丝粒可行，contig唯一对应，但是泛化到其他区间，例如SD可能不行, 已修改
    outfile = open(outfile,'w')
    filter_target_array = {}
    
    
    for i in target_arrays.keys():
        # if i != 'chr10_RagTag_hap1:114448959-115450181':
        #     continue
        sorted_regions = sorted(target_arrays[i],key=lambda x:x[2]) # 按照target asm array的start排序
        # 去重，后一个如果跟前一个重叠，且重叠大于后者的try 30%，去掉后者
        final_pairs = []
        if len(sorted_regions) < 2:
            final_pairs = sorted_regions
            if len(sorted_regions) == 0:
                print(i)
        else:
            final_pairs.append(sorted_regions[0])
            for j in range(len(sorted_regions) - 1):
                if sorted_regions[j + 1][2] > final_pairs[-1][3]:
                    final_pairs.append(sorted_regions[j + 1])
                else:
                    if sorted_regions[j + 1][3] <= final_pairs[-1][3]:
                        continue
                    else:
                        # 检查重叠大小
                        if (final_pairs[-1][3]- sorted_regions[j + 1][2]) / int(sorted_regions[j + 1][1]) > overlap_merge_thr:
                            continue
                        else:
                            final_pairs.append(sorted_regions[j + 1])
        
        # final_pairs, 输出了最终的记录
        contig_table = {} # 建立了一个key:hifiasm contig, value:记录
        for j in final_pairs:
            if j[5] not in contig_table.keys():
                contig_table[j[5]] = [j]
            else:
                contig_table[j[5]].append(j)
  
        
        filter_target_array[i] = {}
        # 遍历contig_table，如果有多个记录保留一个最大的一个，如果该hifiasm contig在该array上是碎裂的，就选择最大的主体
        for j in contig_table.keys():
            record = contig_table[j] # [verkko cen_array,array lenth,array_start,array_end,'+',hifiasm contig name,start,end, hifiasm strand, MQ]
            if len(record) == 1:
                filter_target_array[i][j] = record[0]
            else:
                # 找最长的
                max_length = -1
                max_record = []
                for k in record:
                    if (int(k[7]) - int(k[6])) > max_length:
                        max_length = (int(k[7]) - int(k[6]))
                        max_record = k
                filter_target_array[i][j] = max_record
    
    # 进行过滤操作
    # final_match key:hifiasm-contig, value:array
    final_target_array = {} 
    # 遍历filter_target_array，确定记录的hifiasm-contig是verkko的array之后才加入
    for i in filter_target_array.keys():
        final_target_array[i] = []    
        for j in filter_target_array[i].keys():
            final_target_array[i].append(filter_target_array[i][j])
    
    for i in final_target_array.keys():
        for j in final_target_array[i]:
            record = j
            outfile.write(record[0] + '\t' + record[1] + '\t' + str(record[2]) + '\t' + str(record[3]) + '\t' + record[4] + '\t' + record[5] + '\t' + str(record[6]) + '\t' + str(record[7]) + '\t' +record[8] +'\t' + record[9]+ '\n')
    outfile.close()
 
def buildInitNucflagRegions(matching_records_file,seconda_asm_array_bed,outrepair_table):
    outseconda_asm_array_bed = open(seconda_asm_array_bed,'w')
    outrepair_table = open(outrepair_table,'w')
    
    repair_record = {}
    
    with open(matching_records_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            array = items[0]
            # 不能用contig name作为key
            arrayinfo = array.split(':')
            region = arrayinfo[1].split('-')
            contig = items[5]
            start = int(items[6]) + 1
            end = int(items[7])
            outseconda_asm_array_bed.write(contig + '\t' + str(start) + '\t' + str(end) + '\t' + str(end - start) +  '\n')
            if array not in repair_record.keys():
                repair_record[array] = [[region[0],region[1],'+',contig,str(start),str(end),items[8]]]
            else:
                repair_record[array].append([region[0],region[1],'+',contig,str(start),str(end),items[8]])
                
    outseconda_asm_array_bed.close()
    
    for i in repair_record.keys():
        count = 1
        contig = i.split(':')[0]
        for j in repair_record[i]:
            outrepair_table.write(contig + '\t' + j[0] + '\t' + j[1] + '\t' + j[2] + '\t' + j[3] + '\t' + j[4] + '\t' + j[5] + '\t' + j[6] + '\t' + str(count) + '\n')
            count += 1
    
    outrepair_table.close()

def main():
    # build error array file
    parser = argparse.ArgumentParser(description="Make target region and second assembly region")
    parser.add_argument("-tr", "--target_ref_file", help="", required=True)
    parser.add_argument("-sr", "--second_assembly_file", help="", required=True)
    parser.add_argument("-e", "--error_dir", help="", required=True)
    parser.add_argument("-mt", "--mapping_threads", help="",type=int, default=1)
    parser.add_argument("-mq", "--MAPQ_thr", help="", type=int, default=20)
    parser.add_argument("-bc", "--break_contig_merge_thr", help="", type=int, default=5000)
    parser.add_argument("-ma", "--min_array_thr", help="", type=int, default=100000)
    parser.add_argument("-om", "--overlap_merge_thr", help="", type=float, default=0.3)
    parser.add_argument("-o", "--outdir", help="", required=True)
    
    args = parser.parse_args()
    
    target_ref_file = args.target_ref_file
    second_assembly_file = args.second_assembly_file
    error_dir = args.error_dir
    mapping_threads = args.mapping_threads
    MAPQ_thr = args.MAPQ_thr
    break_contig_merge_thr = args.break_contig_merge_thr
    min_array_thr  = args.min_array_thr
    overlap_merge_thr = args.overlap_merge_thr
    outdir = args.outdir
    
    error_cen_arrays = []
    error_coordinates_file = outdir + '/error_cen_array_list.bed' # output 1
    files = os.listdir(error_dir)
    print(files)
    out_error_coordinates_file = open(error_coordinates_file,'w')
    for i in files:
        if not i.endswith('xls'):
            continue
        array = i.replace('.xls', '')
        error_cen_arrays.append(array)
        info = array.split(':')
        contig = info[0]
        region = info[1].split('-')
        out_error_coordinates_file.write(contig + '\t' + region[0] + '\t' + region[1] + '\n')
    out_error_coordinates_file.close()
        
    
    # 1. get target cen array
    target_cen_array = outdir + '/target_cen_array' # output 2
    if not os.path.exists(target_cen_array):
        os.mkdir(target_cen_array)
    run_getTargetArray(target_ref_file,error_coordinates_file,outdir,target_cen_array)
    
    # 2. mapping cen array to second assembly
    cen_array_fa_file = outdir + '/all.cenarray.fa'
    paf_file = cen_array_fa_file + '.paf'  # output 3
    cmd = 'minimap2 -x asm5 -t ' + str(mapping_threads) +' --eqx --cs --secondary=no -s 25000 -K 8G ' + second_assembly_file + ' ' + cen_array_fa_file + ' > ' + paf_file
    os.system(cmd)
    
    # 3. get matching records from paf
    matching_records_file = outdir + '/cenpairs.xls'  # output 4
    getFilterRecord(paf_file,error_coordinates_file,matching_records_file,MAPQ_thr,break_contig_merge_thr,min_array_thr,overlap_merge_thr)
    
    seconda_asm_array_bed = outdir + '/second_asm_array.nucflag.bed'  # output 5
    repair_table = outdir + '/region_match.all.xls'  # output 6
    buildInitNucflagRegions(matching_records_file,seconda_asm_array_bed,repair_table)
    
    # build second asm array fasta
    second_asm_array = outdir + '/second_asm_array'  # output 7
    if not os.path.exists(second_asm_array):
        os.mkdir(second_asm_array)
    with open(seconda_asm_array_bed,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            cmd = 'samtools faidx ' + second_assembly_file + ' ' + items[0] + ':' + items[1] + '-' + items[2] + ' > ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.tmp.fa'
            os.system(cmd)
            cmd = 'seqkit seq -w 0 ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.tmp.fa' + ' > ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.fa'
            os.system(cmd)
            cmd = 'rm ' + second_asm_array +'/'+  items[0] + ':' + items[1] + '-' + items[2] + '.tmp.fa'
            os.system(cmd)


if __name__ == "__main__":
    main()
    
    