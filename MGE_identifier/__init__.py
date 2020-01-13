
def main():
    from MGE_identifier import identitfy_MGE
    import argparse
    ###Argument handling.
    import argparse
    arg_parser = argparse.ArgumentParser(description='');
    #arg_parser.add_argument("coords_input", help="Directory to show-coords tab-delimited input file.");
    arg_parser.add_argument("-r", "--fasta1", help="reference fasta")
    arg_parser.add_argument("-g", "--gbk", help="reference genbank (to extract annotation of putative MGEs)")
    arg_parser.add_argument("-q", "--fasta2", help="query fasta", nargs='+')
    arg_parser.add_argument("-a", "--algo", help="algorythm to use to compare the genome (megablast, nucmer or promer)", default="nucmer")
    arg_parser.add_argument("-m", "--min_gap_size", help="minimum gap size to consider", default=2000)
    arg_parser.add_argument("-d", "--merge_distance", help="max distance to merge 2 gapped regions", default=2000)
    arg_parser.add_argument("-s", "--samtools_depth", help="samtools depth file", default=False, nargs='+')
    arg_parser.add_argument("-f", "--freq_genomes", help="minimum freq to consider gap position (defaul: 0.9)", default=0.9, type=float)
    args = arg_parser.parse_args()

    test_MGE = identitfy_MGE.MGE(args.fasta1,
                                 args.fasta2,
                                 samtools_depth_files=args.samtools_depth,
                                 run_nucmer=True)
    test_MGE._gap_positions2gap_counts(args.freq_genomes)
    test_MGE._gap_ranges_from_gap_positions()
    test_MGE._merge_close_range(int(args.merge_distance)) # 1000
    test_MGE._filter_small_ranges(int(args.min_gap_size)) # 4000
    if args.gbk:
        test_MGE.extract_annotation(test_MGE.filtered_ranges, args.gbk)
    if args.samtools_depth:
        test_MGE.plot_gap_series_plot(test_MGE.working_dir,
                                      test_MGE.reference_cumulated_length,
                                      test_MGE.mge_table,
                                      args.samtools_depth[0])
    else:
        test_MGE.plot_gap_series_plot(test_MGE.working_dir,
                                      test_MGE.reference_cumulated_length,
                                      test_MGE.mge_table,
                                      False)

    with open("gap_complete_ranges.tab", 'w') as f1:
        for row in test_MGE.range_list:
            data = [str(i) for i in row]
            f1.write('\t'.join(data)+'\n')
    with open("gap_merged_ranges.tab", 'w') as f2:

        for row in test_MGE.merged_ranges:
            data = [str(i) for i in row]
            f2.write('\t'.join(data)+'\n')




