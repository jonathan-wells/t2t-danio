#!/usr/bin/env python3

from collections import defaultdict
from scipy import stats
import pandas as pd
import numpy as np

# TODO: 
# Build TEclass/family dictionary
# Decide on parameters to filter by

class TEData:
    """Collection of relevant information for a TE family extracted from RepeatMasker BED file."""

    def __init__(self, tename, tefam, teclass, chroms=[], starts=[], stops=[], strands=[], percdivs=[]):
        self.tename = tename
        self.tefam = tefam
        self.teclass = teclass
        self.chroms = chroms
        self.starts = starts
        self.stops = stops
        self.percdivs = percdivs
        self.strands = strands

        # These attributes can only be set after RMBed class has instantiated all TEData instances. 
        self.dists = []
        self.tandem = None
        self.enriched_chroms = []
    
    def get_dists(self):
        """Extract distances between fragments using chromosome and start/stop coordinates."""
        assert len(self.starts) == len(self.stops)
        for i in range(1, self.nfragments):
            if self.chroms[i-1] != self.chroms[i]:
                continue
            stop = self.stops[i-1]
            start = self.starts[i]
            self.dists.append(start-stop)
    
    def add_strand(self, strand):
        """Convert +/- strand to boolean values."""
        if strand == '+':
            return 0
        elif strand == '-':
            return 1
        else:
            raise ValueError(f'strand {strand} not recognized.')

    @property
    def nfragments(self):
        """Number of fragments annotated by RepeatMasker."""
        return len(self.starts)

    @property
    def mean_div(self):
        """Mean percent divergence of fragments from consensus sequence."""
        return np.mean(self.percdivs)

    @property
    def median_dist(self):
        """Median distance between fragments. Useful for filtering tandem repeats."""
        if self.dists == []:
            return np.nan
        else:
            return np.median(self.dists)

    @property
    def mean_dist(self):
        """Median distance between fragments. Useful for filtering tandem repeats."""
        if self.dists == []:
            return np.nan
        else:
            return np.mean(self.dists)
    @property
    def median_len(self):
        """Median length of fragments."""
        return np.median([self.stops[i] - self.starts[i] for i in range(self.nfragments)]) 
    
    def assess_chrom_enrichment(self, chrom_counts):
        for chrom in set(self.chroms):
            observed_successes = len([c for c in self.chroms if c == chrom])
            p_value = stats.hypergeom(M=sum(chrom_counts.values()), 
                                      n=chrom_counts[chrom], 
                                      N=len(self.chroms)).sf(observed_successes-1)
            if p_value < 0.01:
                self.enriched_chroms.append(chrom.strip('chr'))

    def assess_tandemicity(self):
        if self.median_dist < 1000:
            self.tandem = True
        else:
            self.tandem = False
    
    def __repr__(self) -> str:
        """Return string containing all available class properties."""
        data = [self.tename,
                self.tefam,
                self.teclass,
                self.nfragments,
                self.mean_div,
                self.median_dist,
                self.mean_dist,
                self.median_len,
                ','.join(self.enriched_chroms),
                self.tandem]
        data = [str(x) for x in data]
        return '\t'.join(data)


class RMBed:
    """Collection of TEData instances compiled from RepeatMasker BED file."""

    def __init__(self, filename):
        self.data = {}
        self.chrom_counts = defaultdict(int)
        
        # Fix broken RepBase labels
        with open('../../data/edta-out/label_dict.txt') as input:
            labels = {l.strip().split('#')[0]: l.strip() for l in input}

        with open(filename) as input:
            for line in input:
                line = line.strip().split('\t')
                chrom, start, stop, tename, swscore, strand, percdiv = line[:7]
                tename = labels.get(tename, tename)
                tename, tefamclass = tename.split('#')
                teclass, tefam = tefamclass.split('/')[0], tefamclass.split('/')[-1]
                
                self.chrom_counts[chrom] += 1

                # Create new TEData instance if TE not yet encountered
                if tename not in self.data.keys():
                    new_te = TEData(tename,
                                    tefam,
                                    teclass,
                                    chroms=[chrom], 
                                    starts=[int(start)],
                                    stops=[int(stop)],
                                    percdivs=[float(percdiv)])
                    new_te.add_strand(strand)
                    self.data[tename] = new_te
                
                else: # Add to existing
                    self.data[tename].chroms.append(chrom)
                    self.data[tename].starts.append(int(start))
                    self.data[tename].stops.append(int(stop))
                    self.data[tename].add_strand(strand)
                    self.data[tename].percdivs.append(float(percdiv))

        # Sort TEs by name and calculate secondary attributes, e.g. chromosome enrichment
        self._sorted_data = sorted(self.data.values(), key=lambda x: x.tename)
        for tedata in self._sorted_data:
            tedata.get_dists()
            tedata.assess_chrom_enrichment(self.chrom_counts)
            tedata.assess_tandemicity()

    def to_tsv(self, filename):
        """Write summary of RMBed data to file."""
        with open(filename, 'w') as output:
            output.write('tename\ttefam\tteclass\tnfragments\tmean_div\tmedian_dist\tmean_dist\tmedian_len\tchrom_enriched\tistandem\n')
            for tedata in self:
                # str(tedata) is doing a lot of work here: see TEData.__repr__
                line = str(tedata)
                output.write(f'{line}\n')


    def __getitem__(self, index):
        return self._sorted_data[index]

    def __len__(self):
        return len(self.data.keys())

    def __iter__(self):
        for tedata in self._sorted_data:
            yield tedata



def filter_summary_file(filename):
    rm_df = pd.read_csv(filename, sep='\t')
    print(rm_df.head(10))


if __name__ == "__main__":
    # rmbed = RMBed('../../data/edta-out/fish11t2t.fasta.mod.out.bed')
    # rmbed.to_tsv('../../data/edta-out/parsed_rmbed.txt')
    filter_summary_file('../../data/edta-out/parsed_rmbed.txt')
