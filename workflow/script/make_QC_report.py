import os
import argparse

def main():

    parser = argparse.ArgumentParser(description="Compare errors before and after repair")
    parser.add_argument("-r", "--region_file", help="", required=True)
    parser.add_argument("-tn", "--target_nucflag", help="", required=True)
    parser.add_argument("-an", "--final_nucflag", help="", required=True)
    parser.add_argument("-o", "--qc_report", help="", required=True)
    
    args = parser.parse_args()
    
    region_file = args.region_file
    target_nucflag = args.target_nucflag
    final_nucflag = args.final_nucflag
    qc_report = args.qc_report

    repaired_region_pairs = {}
    with open(region_file,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            items = line.split('\t')
            new_region = items[0] + ':' + items[1] + '-' + items[2]
            old_region = items[0] + ':' + items[4] + '-' + items[5]
            repaired_region_pairs[new_region] = old_region

    outfile = open(qc_report,'w')
    outfile.write('Region\tInit_error_base_num\tFinal_error_base_num\tError_reduction\n')
    all_target_error_num = 0
    all_final_error_num = 0
    for i in repaired_region_pairs.keys():
        new_region = i
        old_region = repaired_region_pairs[i]
        target_file = target_nucflag + '/' + old_region + '.xls'
        target_error_num = 0
        with open(target_file,'r') as f:
            while True:
                line = f.readline()[:-1]
                if not line:
                    break
                items = line.split('\t')
                target_error_num += (int(items[1]) - int(items[0]))
        all_target_error_num += target_error_num
        final_file = final_nucflag + '/' + new_region + '.xls'
        final_error_num = 0
        if not os.path.exists(final_file):
            final_error_num = 0
        else:
            with open(final_file,'r') as f:
                while True:
                    line = f.readline()[:-1]
                    if not line:
                        break
                    items = line.split('\t')
                    
                    final_error_num += (int(items[1]) - int(items[0]))
        all_final_error_num += final_error_num
        error_reduction = target_error_num - final_error_num
        outfile.write(i.split('.xls')[0] + '\t' + str(target_error_num) + '\t' + str(final_error_num) + '\t' + str(error_reduction) + '\n')
    outfile.write('All\t' + str(all_target_error_num) + '\t' + str(all_final_error_num) + '\t' + str(all_target_error_num - all_final_error_num) + '\n')
    outfile.close()
    
    
    
    

if __name__ == "__main__":
    main()
    
    
    
    
    