# -*- coding: utf-8 -*-
# v.140216

# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the docs/COPYING file.

import sys
import numpy
import hashlib


class FastaOutput:
    def __init__(self, output_file_path):
        self.output_file_path = output_file_path
        self.output_file_obj = open(output_file_path, 'w')

    def store(self, entry, split=True, store_frequencies=True):
        if entry.unique and store_frequencies:
            self.write_id('%s|%s' % (entry.id, 'frequency:%d' % len(entry.ids)))
        else:
            self.write_id(entry.id)

        self.write_seq(entry.seq, split)

    def write_id(self, id_):
        self.output_file_obj.write('>%s\n' % id_)

    def write_seq(self, seq, split=True):
        seq = seq.decode() if isinstance(seq, bytes) else seq
        if split:
            seq = self.split(seq)
        self.output_file_obj.write('%s\n' % seq)

    def split(self, seq, piece_length=80):
        seq = seq.decode() if isinstance(seq, bytes) else seq

        ticks = list(range(0, len(seq), piece_length)) + [len(seq)]
        return '\n'.join([seq[ticks[x]:ticks[x + 1]] for x in range(0, len(ticks) - 1)])

    def close(self):
        self.output_file_obj.close()


class ReadFasta:
    def __init__(self, f_name):
        self.ids = []
        self.sequences = []

        self.fasta = SequenceSource(f_name)

        while next(self.fasta):
            if self.fasta.pos % 1000 == 0 or self.fasta.pos == 1:
                sys.stderr.write('\r[fastalib] Reading FASTA into memory: %s' % self.fasta.pos)
                sys.stderr.flush()
            self.ids.append(self.fasta.id)
            self.sequences.append(self.fasta.seq)

        sys.stderr.write('\n')

    def close(self):
        self.fasta.close()


class SequenceSource(object):
    def __init__(self, fasta_file_path, lazy_init=True, unique=False, allow_mixed_case=False, progress_obj=None):
        self.fasta_file_path = fasta_file_path
        self.name = None
        self.lazy_init = lazy_init
        self.allow_mixed_case = allow_mixed_case
        self.progress_obj=progress_obj

        self.pos = 0
        self.id = None
        self.md5id = None
        self.seq = None
        self.ids = []

        self.unique = unique
        self.total_unique = 0
        self.max_freq = 0
        self.unique_hash_dict = {}
        self.unique_hash_list = []
        self.unique_next_hash = 0

        self.file_pointer = open(self.fasta_file_path)
        self.file_pointer.seek(0)

        if self.lazy_init:
            self.number_seq = None
        else:
            self.number_seq = len([1 for ln in self.file_pointer.readlines() if ln.startswith('>')])
            self.reset()

        if self.unique:
            self.init_unique_hash()

    def init_unique_hash(self):
        if self.progress_obj:
            self.progress_obj.append('Hashing the reads into memory: %s' % self.pos)
            
        while self.next_regular():
            if self.progress_obj and self.pos % 5000 == 0:
                self.progress_obj.update('Hashing the reads into memory: %s' % self.pos)

            seq_hash = hashlib.sha1(self.seq.upper()).hexdigest()
            if seq_hash in self.unique_hash_dict:
                self.unique_hash_dict[seq_hash]['ids'].append(self.id)
                self.unique_hash_dict[seq_hash]['count'] += 1
            else:
                self.unique_hash_dict[seq_hash] = {'id': self.id,
                                                   'ids': [self.id],
                                                   'seq': self.seq,
                                                   'count': 1}

        sorted_descending_by_counts = sorted([(self.unique_hash_dict[seq_hash]['count'], seq_hash)
                                              for seq_hash in self.unique_hash_dict], reverse=True)
        self.unique_hash_list = [hs for _, hs in sorted_descending_by_counts]
        self.max_freq = sorted_descending_by_counts[0][0]
        self.total_unique = len(self.unique_hash_dict)
        self.reset()

    def __iter__(self):
        return self

    def __next__(self):
        if self.unique:
            return self.next_unique()
        else:
            return self.next_regular()

    def next_unique(self):
        if self.unique:
            if self.total_unique > 0 and self.pos < self.total_unique:
                self.md5id = self.unique_hash_list[self.pos]
                hash_entry = self.unique_hash_dict[self.md5id]

                self.pos += 1
                self.seq = hash_entry['seq'] if self.allow_mixed_case else hash_entry['seq'].upper()
                self.id = hash_entry['id']
                self.ids = hash_entry['ids']

                return True
            else:
                return False
        else:
            return False

    def next_regular(self):
        self.seq = None
        self.id = self.file_pointer.readline()[1:].strip()
        sequence = ''

        while 1:
            line = self.file_pointer.readline()
            if not line:
                if len(sequence):
                    self.seq = sequence.encode("utf-8")
                    self.pos += 1
                    return True
                else:
                    return False
            if line.startswith('>'):
                self.file_pointer.seek(self.file_pointer.tell() - len(line))
                break
            sequence += line.strip()

        self.seq = sequence if self.allow_mixed_case else sequence.upper()
        self.seq = self.seq.encode("utf-8")
        self.pos += 1
        return True

    def get_seq_by_read_id(self, read_id):
        self.reset()
        while next(self):
            if self.id == read_id:
                return self.seq

        return False

    def close(self):
        self.file_pointer.close()

    def reset(self):
        self.pos = 0
        self.id = None
        self.seq = None
        self.ids = []
        self.file_pointer.seek(0)

    def visualize_sequence_length_distribution(self, title, dest=None, max_seq_len=None,
                                               xtickstep=None, ytickstep=None):
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        sequence_lengths = []

        self.reset()

        while next(self):
            if self.pos % 10000 == 0 or self.pos == 1:
                sys.stderr.write('\r[fastalib] Reading: %s' % self.pos)
                sys.stderr.flush()
            sequence_lengths.append(len(self.seq))

        self.reset()

        sys.stderr.write('\n')

        if not max_seq_len:
            max_seq_len = max(sequence_lengths) + (int(max(sequence_lengths) / 100.0) or 10)

        seq_len_distribution = [0] * (max_seq_len + 1)

        for ln in sequence_lengths:
            seq_len_distribution[ln] += 1

        fig = plt.figure(figsize=(16, 12))
        plt.rcParams.update({'axes.linewidth': 0.9})
        plt.rc('grid', color='0.50', linestyle='-', linewidth=0.1)

        gs = gridspec.GridSpec(10, 1)

        ax1 = plt.subplot(gs[0:8])
        plt.grid(True)
        plt.subplots_adjust(left=0.05, bottom=0.03, top=0.95, right=0.98)

        plt.plot(seq_len_distribution, color='black', alpha=0.3)
        plt.fill_between(list(range(0, max_seq_len + 1)), seq_len_distribution, y2=0, color='black', alpha=0.15)
        plt.ylabel('number of sequences')
        plt.xlabel('sequence length')

        if xtickstep is None:
            xtickstep = (max_seq_len / 50) or 1

        if ytickstep is None:
            ytickstep = max(seq_len_distribution) / 20 or 1

        plt.xticks(list(range(xtickstep, max_seq_len + 1, xtickstep)), rotation=90, size='xx-small')
        plt.yticks(list(range(0, max(seq_len_distribution) + 1, ytickstep)),
                   [y for y in range(0, max(seq_len_distribution) + 1, ytickstep)],
                   size='xx-small')
        plt.xlim(xmin=0, xmax=max_seq_len)
        plt.ylim(ymin=0, ymax=max(seq_len_distribution) + (max(seq_len_distribution) / 20.0))

        plt.figtext(0.5, 0.96, '%s' % title, weight='black', size='xx-large', ha='center')

        ax1 = plt.subplot(gs[9])
        plt.grid(False)
        plt.yticks([])
        plt.xticks([])
        plt.text(0.02, 0.5, 'total: %s / mean: %.2f / std: %.2f / min: %s / max: %s'
                 % (len(sequence_lengths),
                    numpy.mean(sequence_lengths), numpy.std(sequence_lengths),
                    min(sequence_lengths),
                    max(sequence_lengths)),
                 va='center', alpha=0.8, size='x-large')

        if dest is None:
            dest = self.fasta_file_path

        try:
            plt.savefig(dest + '.pdf')
        except:
            plt.savefig(dest + '.png')

        try:
            plt.show()
        except:
            pass

        return


class QualSource:
    def __init__(self, quals_file_path, lazy_init=True):
        self.quals_file_path = quals_file_path
        self.name = None
        self.lazy_init = lazy_init

        self.pos = 0
        self.id = None
        self.quals = None
        self.quals_int = None
        self.ids = []

        self.file_pointer = open(self.quals_file_path)
        self.file_pointer.seek(0)

        if self.lazy_init:
            self.total_quals = None
        else:
            self.total_quals = len([l for l in self.file_pointer.readlines() if l.startswith('>')])
            self.reset()

    def __next__(self):
        self.id = self.file_pointer.readline()[1:].strip()
        self.quals = None
        self.quals_int = None

        qualscores = ''

        while 1:
            line = self.file_pointer.readline()
            if not line:
                if len(qualscores):
                    self.quals = qualscores.strip()
                    self.quals_int = [int(q) for q in self.quals.split()]
                    self.pos += 1
                    return True
                else:
                    return False
            if line.startswith('>'):
                self.file_pointer.seek(self.file_pointer.tell() - len(line))
                break
            qualscores += ' ' + line.strip()

        self.quals = qualscores.strip()
        self.quals_int = [int(q) for q in self.quals.split()]
        self.pos += 1

        return True

    def close(self):
        self.file_pointer.close()

    def reset(self):
        self.pos = 0
        self.id = None
        self.quals = None
        self.quals_int = None
        self.ids = []
        self.file_pointer.seek(0)


if __name__ == '__main__':
    fasta = SequenceSource(sys.argv[1])
    fasta.visualize_sequence_length_distribution(title=sys.argv[2] if len(sys.argv) == 3 else 'None')
