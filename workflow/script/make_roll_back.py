import argparse
import os


def main():
    #
    # 检查对应区间原来的错误总碱基和当前的错误总碱基，如果没有减少则回退
    #
    parser = argparse.ArgumentParser(description="Compare errors before and after repair")
    parser.add_argument("-r", "--repair_region_file", help="", required=True)
    parser.add_argument("-or", "--old_region_error_dir", help="", required=True)
    parser.add_argument("-nr", "--new_region_error_dir", help="", required=True)
    parser.add_argument("-o", "--outfile", help="", required=True)
    
    args = parser.parse_args()
    
    repair_region_file = args.repair_region_file
    old_region_error_dir = args.old_region_error_dir
    new_region_error_dir = args.new_region_error_dir
    outfile = args.outfile
    
    # repair_region_file = '/project/logsdon_shared/projects/HGSVC3/AssemblyRepairer_new/Snakemake/results_target/repaired.bed'
    # old_region_error_dir = '/project/logsdon_shared/projects/HGSVC3/AssemblyRepairer_new/Snakemake/results_target/nucflag_init_target/merge_errors'
    # new_region_error_dir = '/project/logsdon_shared/projects/HGSVC3/AssemblyRepairer_new/Snakemake/results_target/nucflag_repaired_asm/merge_errors'
    # outdir = '/project/logsdon_shared/projects/HGSVC3/AssemblyRepairer_new/Snakemake/results_target'
    
    
    repaired_regions = {}
    repaired_state = {}
    with open(repair_region_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            new_region = items[0] + ':' + items[1] + '-' + items[2]
            old_region = items[0] + ':' + items[4] + '-' + items[5]
            repaired_regions[new_region] = old_region
    
    error_base_num = {}
    
    for i in repaired_regions.keys():
        old_file = old_region_error_dir + '/' + repaired_regions[i] + '.xls'
        new_file = new_region_error_dir + '/' + i + '.xls'
        old_error_num = 0
        new_error_num = 0
        if os.path.exists(old_file):
            with open(old_file,'r') as f:
                while True:
                    line = f.readline()[:-1]
                    if not line:
                        break
                    items = line.split('\t')
                    old_error_num += (int(items[1]) - int(items[0]))
        if os.path.exists(new_file):
            with open(new_file,'r') as f:
                while True:
                    line = f.readline()[:-1]
                    if not line:
                        break
                    items = line.split('\t')
                    new_error_num += (int(items[1]) - int(items[0]))
        error_base_num[i] = [old_error_num,new_error_num]
        
    for i in error_base_num.keys():
        if error_base_num[i][0] <= error_base_num[i][1]:
            repaired_state[i] = ['error',repaired_regions[i],error_base_num[i][0],error_base_num[i][1]]
        else:
            repaired_state[i] = ['success',repaired_regions[i],error_base_num[i][0],error_base_num[i][1]]
        
    outRollback_file_name = outfile
    outRollback_file = open(outRollback_file_name,'w')
    for i in repaired_state.keys():
        outRollback_file.write(i + '\t' + repaired_state[i][1] + '\t' + repaired_state[i][0] + '\t' + str(repaired_state[i][2]) + '\t' + str(repaired_state[i][3]) + '\n' )
        
    outRollback_file.close()        
                
    
    


if __name__ == '__main__':
    main()