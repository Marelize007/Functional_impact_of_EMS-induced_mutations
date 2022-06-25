# The follow python script can be used to filter out EMS induced mutations. The input file is a VCF file containing the tentative mutations of all treatment lines per isolate.
# The output files are a VCF file containing the EMS induced mutations that passed the filtering criteria, a text file with a break down of the types of base substitutions, as well as a text file containing all the mutations that did not pass the filtering criteria.

file_path = "input.vcf"
output_path = "output.vcf"
tstv_summary_path = "tstv_ratio.txt"
failed_path = "failed.txt"

min_count = 2     # number of forward and revers reads required to make a mutation call
min_sample_consensus = 12  # minimun amount of samples needed to produce a consensus

sample_column_list = [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]  # enter all column index of samples
column_of_interest = 9 # sample you are detecting mutations for

progress = 0
# make dict to store counts of all base substitution
base_sub_dict = dict(zip(['AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TG'], [0]*12))
with open(failed_path, 'w+') as failed_file, open(output_path, 'w+') as output_file, open(file_path, 'r') as input_file:
            for line in input_file:

                # let's you know the code did not freeze
                progress += 1
                if progress % 1000 == 0:
                    print(progress)

                # prints out header to output
                if line.startswith("#"):
                    output_file.write(line)
                else:
                    info_pos = []
                    genotype = []
                    write_line = True

                    columns = line.strip().split()

                    # stores all possible alleles; index will correspond with genotype call
                    ref_alt_alleles = [columns[3]]  # stores the ref allele
                    ref_alt_alleles.extend(columns[4].split(','))  # stores all alt alleles

                    # stores "[0/1, 36,22, 22,12, 14,10, 58, 255,0,255, 1, 1]" in info_pos
                    for sample_column in sample_column_list:
                        info_pos.append(columns[sample_column].split(":"))

                    if len(info_pos[0]) > 2:  
                        for sample in info_pos:

                            fwd_rev_pos = []
                            # stores fwd rev reads as list of ints [rf, rr, af, ar]
                            fwd_list, rev_list = [list(map(int, x.split(','))) for x in sample[4:6]]
                            for fwd, rev in zip(fwd_list, rev_list):
                                fwd_rev_pos.extend([fwd, rev])
                            # vets each read count against marelize's standards
                            good_enough_for_marelize = [read_count >= min_count for read_count in fwd_rev_pos]

                            fwd_rev_pair_list = []
                            for i in range(int(len(good_enough_for_marelize) / 2)):
                                fwd_rev_pair_list.append([i * 2, i * 2 + 1])
                                
                            if sum(good_enough_for_marelize) == 2:
                                if sum([fwd_rev_pos[i] for i, x in enumerate(good_enough_for_marelize) if not x]) != 0:
                                    write_line = False
                                if [idx for idx, good_enough in enumerate(good_enough_for_marelize) if good_enough] in fwd_rev_pair_list:
                                    homozygous_allele = str(fwd_rev_pair_list.index([idx for idx, good_enough in enumerate(good_enough_for_marelize) if good_enough]))
                                    genotype.append(homozygous_allele + "/" + homozygous_allele)
                                else:
                                    genotype.append("./.")
                            elif sum(good_enough_for_marelize) == 4:
                                idx_good_enough = [idx for idx, good_enough in enumerate(good_enough_for_marelize) if good_enough]
                                idx_good_enough = [idx_good_enough[:2], idx_good_enough[2:]]
                                heterozygous_allele = []
                                for pair in idx_good_enough:
                                    if pair in fwd_rev_pair_list:
                                        heterozygous_allele.append(str(fwd_rev_pair_list.index(pair)))
                                if len(heterozygous_allele) == 2:
                                    genotype.append(heterozygous_allele[0] + "/" + heterozygous_allele[1])
                                else:
                                    genotype.append("./.")
                            else:
                                genotype.append("./.")

                    # removes all ambigous genotypes and stores it in called_gt
                    called_gt = [gt for gt in genotype if gt != "./."]
                    consensus_gt = ""
                    if called_gt:
                        # finds the consensus genotype and checks if it is within the relaxation range
                        consensus_gt = max(set(called_gt), key=called_gt.count)

                        if genotype.count(consensus_gt) < min_sample_consensus:
                            write_line = False

                        if consensus_gt[0] != consensus_gt[-1]:
                            write_line = False

                    else:
                        write_line = False


                    # finds index of column of interest in sample_column_list
                    coi_idx = sample_column_list.index(column_of_interest)

                    # checks that the genotype of the sample of interest is unique
                    if write_line:
                        coi_genotype = genotype.pop(coi_idx)
                        called_gt = [gt for gt in genotype if gt != "./."]
                        # makes sure coi_genotype is different from consensus
                        if coi_genotype == consensus_gt:
                            write_line = False
                            
                        else:
                            # only use mutations that change from homozygous to heterzygous
                            if coi_genotype[0] != coi_genotype[-1] and consensus_gt[0] == consensus_gt[-1]:
                                # records the type of base substitution
                                coi_bases = coi_genotype.split('/')
                                consen_bases = consensus_gt.split('/')
                                for base in coi_bases:
                                    try:
                                        consen_bases.remove(base)
                                    except ValueError:
                                        base_sub_type = ref_alt_alleles[int(consen_bases[0])] + ref_alt_alleles[int(base)]
                                        base_sub_dict[base_sub_type] +=  1
                            else:
                                write_line = False

                            # reinserts coi gt into genotype for output later
                            genotype.insert(coi_idx, coi_genotype)

                    if write_line:
                        # replaces the old genotype call with our own genotype call and reassembly the line to output
                        for sample_idx, new_genotype in enumerate(genotype):
                            info_pos[sample_idx][0] = new_genotype
                        for sample_column, sample in zip(sample_column_list, info_pos):
                            columns[sample_column] = ':'.join(sample)
                        line = '\t'.join(columns)

                        output_file.write(line + '\n')
                    else:
                        # we do not attempt to call genotype for lines that fail, so it will have the old genotype call
                        failed_file.write(line + '\n')

with open(tstv_summary_path, 'w+') as tstv_output:
    ts_all = ['AG', 'GA', 'CT', 'TC']
    tv_all = ['AC', 'AT', 'CA', 'CG', 'GC', 'GT', 'TA', 'TG']
    ts_counts = sum([base_sub_dict[ts] for ts in ts_all])
    tv_counts = sum([base_sub_dict[tv] for tv in tv_all])
    try:
        tstv_ratio = ts_counts / tv_counts
        tstv_output.write("ts/tv ratio:\t" + "%.2f" % tstv_ratio + "\n-----\n")
    except ZeroDivisionError:
        tstv_output.write("ts/tv ratio:\t" + "NA" + "\n-----\n")
    for base_sub_type, count in base_sub_dict.items():
        tstv_output.write(base_sub_type[0] + '>' + base_sub_type[1] + ':\t' + str(count) + '\n')                         
                              
