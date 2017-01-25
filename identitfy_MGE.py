#!/usr/bin/python

import nucmer_utility

class MGE:
    def __init__(self, reference_genome, query_genomes, samtools_depth_files=False):

        import numpy as np
        import pandas as pd
        from Bio import SeqIO

        if samtools_depth_files:
            self.samtools_depth_files = [pd.read_csv(i, sep="\t",names=['chr','position','depth']) for i in samtools_depth_files]
        else:
            self.samtools_depth_files = False


        self.reference_genome = reference_genome
        self.reference_genome_records = [i for i in SeqIO.parse(open(self.reference_genome, 'r'), 'fasta')]
        self.reference_cumulated_length = self._get_cumulated_fasta_length(self.reference_genome_records)
        self.cumulated_sequences = self._get_cumulated_fasta_seq(self.reference_genome_records)
        self.filtered_ranges = []


        if not isinstance(query_genomes, list):
            self.query_genomes = list(query_genomes)
        else:
            self.query_genomes = query_genomes

        self.n_query_genomes = len(self.query_genomes)

        self.query_genome2gap_pistions =  {}

        self.coord_files, self.delta_files = nucmer_utility.execute_promer(self.reference_genome,
                                                            self.query_genomes,
                                                            coords=False)
        self.contig_add = nucmer_utility.get_contigs_coords(self.reference_genome)


        for delta_file in self.delta_files:

            print 'contig adds', self.contig_add

            contig2start_stop_list = nucmer_utility.delta_file2start_stop_list(delta_file,
                                       contigs_add=self.contig_add,
                                       algo='nucmer',
                                       minimum_identity=85)
            gap_data = nucmer_utility.get_gaps_from_start_stop_lists(contig2start_stop_list,
                                           contigs_add=self.contig_add,
                                           min_gap_size=1000)

            query_genome = delta_file.split('.')[0]
            start_list = []
            stop_list = []
            for contig in gap_data:
                for gap in gap_data[contig]:
                    start_list.append(gap[0])
                    stop_list.append(gap[1])
                data = pd.DataFrame({'start': start_list, 'stop': stop_list })
                data.start = data.start.astype(np.int64)
                data.stop = data.stop.astype(np.int64)
                data_sort = data.sort(columns=["start"])
                data_sort.to_csv("%s_gaps.csv" % delta_file.split('.')[0])

                self.query_genome2gap_pistions[query_genome] = data_sort

    def _gap_positions2gap_counts(self, cutoff=0.9):
        '''
        counts number of gaps for each position in the genome
        calculate fraction of gaps out of the number of query genomes
        :return:
        '''



        self.genome_counts = [0]*self.reference_cumulated_length
        for query_genome in self.query_genome2gap_pistions:
            gaps = self.query_genome2gap_pistions[query_genome]
            for i in range(0, len(gaps['start'])-1):
                for gap_position in range(gaps['start'][i],gaps["stop"][i]):
                    self.genome_counts[gap_position] +=1

        self.genome_counts_fractions = [i/float(self.n_query_genomes) for i in self.genome_counts]
        self.gap_positions = [n for n,i in enumerate(self.genome_counts_fractions) if i>=cutoff]

        print 'length gap positions', len(self.gap_positions)

    def _gap_ranges_from_gap_positions(self):
        self.range_list = []
        new_range = [self.gap_positions[0],self.gap_positions[0]]
        print self.gap_positions[0:5000]
        for position in range(1,len(self.gap_positions)):
            if self.gap_positions[position] == self.gap_positions[position-1] + 1:
                new_range[1]+=1
            else:
                self.range_list.append(new_range)
                new_range = [self.gap_positions[position], self.gap_positions[position]]
        print len(self.range_list)
        print self.range_list

    def _merge_close_range(self, max_distance):
        self.merged_ranges = []
        previous_range = self.range_list[0]
        for one_range in range(1,len(self.range_list)):
            start_range = self.range_list[one_range][0]
            stop_range = self.range_list[one_range][1]

            start_previous_range = previous_range[0]
            stop_previous_range = previous_range[1]

            # if the 2 ranges are closer than <max distance>, merge the 2 range
            # and redefine <previous grange>
            if start_range-stop_previous_range <= max_distance:
                if one_range < len(self.range_list)-1:
                    previous_range = [start_previous_range, stop_range]
                else:
                    # last range, add the merged range
                    previous_range = [start_previous_range, stop_range]
                    #print 'adding:', previous_range
                    self.merged_ranges.append(previous_range)
            else:
                self.merged_ranges.append(previous_range)
                if one_range < len(self.range_list)-1:
                    # kepp the previous range in memory and set a new one
                    previous_range = self.range_list[one_range]
                else:
                    # last range, add the current range anyway
                    self.merged_ranges.append(self.range_list[one_range])


    def extract_annotation(self, ranges, genbank_file):
        from Bio import SeqIO

        with open(genbank_file, 'r') as g:
            records = [i for i in SeqIO.parse(g, 'genbank')]
            if len(records) > 1:
                merged_record = ''
                for record in records:
                    merged_record+=record
                records = merged_record
            else:
                records =records[0]
        for n, mge_range in enumerate(ranges):
            with open('MGE_%s.gbk' % n, 'w') as o:
                mge_record = records[mge_range[0]:mge_range[1]]
                mge_record.id = 'MGE_%s' % n
                mge_record.accession = 'MGE_%s' % n
                mge_record.name = 'MGE_%s' % n
                print 'mge %s, length: %s, n features: %s' % (n,
                                                              len(mge_record.seq),
                                                              len(mge_record.features))
                SeqIO.write(mge_record,o, 'genbank')



    def _filter_small_ranges(self, minimal_size=4000):
        from Bio.SeqUtils import GC
        from Bio import Seq
        import statistics
        import re
        i=1
        print len(self.merged_ranges)

        for range in self.merged_ranges:
            start = range[0]
            stop = range[1]
            seq = self.cumulated_sequences[start:stop].seq
            if 'N' in str(seq):
                gaps = 'YES'
                # remove 'nnn' before calculating GC
                seq_bis = Seq.Seq(re.sub('N', '', str(seq)))
                gc_content = round(GC(seq_bis),2)
            else:
                gaps = 'NO'
                gc_content = round(GC(seq), 2)
            if self.samtools_depth_files is not False:
                median_depth = []
                sd_depth = []
                for depth_file in self.samtools_depth_files:
                    median_depth.append(statistics.median(depth_file.loc[start:stop]['depth']))
                    sd_depth.append(round(statistics.stdev(depth_file.loc[start:stop]['depth']),2))


            if range[1]-range[0] >= minimal_size:
                self.filtered_ranges.append(range)
                mge_data = 'MGE_%s\t%s\t%s\t%s\t%s\t%s' % (i,
                                                              start,
                                                              stop,
                                                              range[1]-range[0],
                                                              gaps,
                                                              gc_content)
                if self.samtools_depth_files:
                    for n, depth in enumerate(self.samtools_depth_files):
                        mge_data+="\t%s\t%s" % (median_depth[n], sd_depth[n])
                print mge_data
            # MGE counts
                i+=1

        n_contigs = len(str(self.cumulated_sequences.seq).split('N'*200))
        gappless_chromosome = Seq.Seq(re.sub('N', '', str(self.cumulated_sequences.seq)))
        gc_content_chromosme = round(GC(gappless_chromosome),2)
        if self.samtools_depth_files is not False:
            median_depth = []
            sd_depth = []
            for depth_file in self.samtools_depth_files:
                median_depth.append(statistics.median(depth_file.loc[:]['depth']))
                sd_depth.append(round(statistics.stdev(depth_file.loc[:]['depth']),2))


        chr_data = 'CHROMOSOME\t%s\t%s\t%s\t%s\t%s' % (0,
                                                      len(gappless_chromosome),
                                                      len(gappless_chromosome),
                                                      n_contigs,
                                                      gc_content_chromosme)
        if self.samtools_depth_files:
            for n, depth in enumerate(self.samtools_depth_files):
                chr_data+="\t%s\t%s" % (median_depth[n], sd_depth[n])
        print chr_data


    def _get_cumulated_fasta_length(self, records):
        length = 0
        for record in records:
            length+=len(record.seq)
        return length

    def _get_cumulated_fasta_seq(self, records):
        merged_record = ''
        for record in records:
            merged_record+=record
        return merged_record

if __name__ == '__main__':
    ###Argument handling.
    import argparse
    arg_parser = argparse.ArgumentParser(description='');
    #arg_parser.add_argument("coords_input", help="Directory to show-coords tab-delimited input file.");
    arg_parser.add_argument("-r", "--fasta1", help="reference fasta")
    arg_parser.add_argument("-g", "--gbk", help="reference genbank (to extract annotation of putative MGEs)")
    arg_parser.add_argument("-q", "--fasta2", help="query fasta", nargs='+')
    arg_parser.add_argument("-a", "--algo", help="algorythm to use to compare the genome (megablast, nucmer or promer)", default="nucmer")
    arg_parser.add_argument("-m", "--min_gap_size", help="minimum gap size to consider", default=1000)
    arg_parser.add_argument("-s", "--samtools_depth", help="samtools depth file", default=False, nargs='+')
    arg_parser.add_argument("-f", "--freq_genomes", help="minimum freq to consider gap position (defaul: 0.9)", default=0.9, type=float)
    args = arg_parser.parse_args()
    print args.freq_genomes, type(args.freq_genomes)
    test_MGE = MGE(args.fasta1, args.fasta2, samtools_depth_files=args.samtools_depth)
    test_MGE._gap_positions2gap_counts(args.freq_genomes)
    test_MGE._gap_ranges_from_gap_positions()
    test_MGE._merge_close_range(1000) # 10000
    test_MGE._filter_small_ranges(4000)
    if args.gbk:
        test_MGE.extract_annotation(test_MGE.filtered_ranges, args.gbk)


