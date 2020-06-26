# -*- coding: utf-8

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
import time  # agc
import gc
import os
import numpy
import operator

from Oligotyping.lib import fastalib as u
from Oligotyping.lib.entropy import entropy
from Oligotyping.utils.utils import ConfigError


class Topology:
    def __init__(self, nodes_output_directory=None):
        self.nodes = {}
        self.alive_nodes = []
        self.zombie_nodes = []
        self.final_nodes = []

        # this is a temporary register to store node ids that needs to be
        # checked by certain functions
        self.standby_bin = []

        self.nodes_output_directory = nodes_output_directory

        self.outliers = {}
        self.outlier_reasons = []
        self.next_available_node_id = None

        self.frequency_of_the_most_abundant_read = None
        self.average_read_length = None
        self.alignment_length = None

        self.logger = None

    def get_new_node_id(self):
        if not self.next_available_node_id:
            new_node_id = 1
            self.next_available_node_id = new_node_id + 1
        else:
            new_node_id = self.next_available_node_id
            self.next_available_node_id += 1

        return '%.9d' % new_node_id

    def add_new_node(self, node_id, unique_read_objects_list, root=False, parent_id=None):
        if not self.nodes_output_directory:
            raise ConfigError("Nodes output directory has to be declared before adding new nodes")

        # node = Node(node_id, self.nodes_output_directory)
        node = Node(node_id, self.nodes_output_directory, self.logger)  # agc
        node.reads = unique_read_objects_list
        node.size = sum([read.frequency for read in node.reads])
        node.pretty_id = self.get_pretty_id(node_id)

        if parent_id:
            parent = self.nodes[parent_id]
            parent.children.append(node.node_id)
            node.level = parent.level + 1
            node.parent = parent.node_id

        node.refresh()

        if root:
            # things to initialize if this is the root node
            node.level = 0

            self.frequency_of_the_most_abundant_read = max([read.frequency for read in node.reads])
            self.alignment_length = len(node.reads[0].seq)

            # store the average read length. this is a terrible approximation,
            # but it is way better to go through all reads which are most probably
            # have the same length anyway (or minus/plus 2-3 nt at worst).
            self.average_read_length = len(node.reads[0].seq.replace(b'-', b''))

        self.nodes[node_id] = node

        return node

    def get_pretty_id(self, node_id):
        pretty_id = node_id
        while 1:
            if pretty_id[0] == '0':
                pretty_id = pretty_id[1:]
            else:
                break
        return pretty_id

    def get_node(self, node_id):
        return self.nodes[node_id]

    def print_node(self, node_id):
        node = self.nodes[node_id]
        print()
        print('Node "%s"' % node)
        print('---------------------------------')
        print('Alive     : %s' % (not node.killed))
        print('Dirty     : %s' % node.dirty)
        print('Size      : %d' % node.size)
        print('Parent    : %s' % node.parent)
        print('Children  :', node.children)
        print()

    def get_final_count(self):
        return sum([sum([read.frequency for read in self.nodes[node_id].reads]) for node_id in self.final_nodes])

    def get_siblings(self, node_id):
        node = self.nodes[node_id]
        siblings = []

        parent_nodes_to_analyze = [self.get_node(node.parent)]

        while len(parent_nodes_to_analyze):
            parent = parent_nodes_to_analyze.pop(0)

            for candidate in [self.get_node(c) for c in parent.children if c != node.node_id]:
                if candidate.children:
                    parent_nodes_to_analyze.append(candidate)
                else:
                    siblings.append((candidate.size, candidate.node_id))

        # sort siblings least abundant to most
        siblings.sort()

        return [s[1] for s in siblings]

    def remove_node(self, node_id, store_content_in_outliers_dict=False, reason=None):
        node = self.nodes[node_id]
        parent = self.nodes[node.parent]

        self.logger.info('topology: remove node "%s" (of parent "%s")' % (node_id, node.parent))

        parent.children.remove(node_id)
        parent.size -= node.size

        if store_content_in_outliers_dict:
            for read in node.reads:
                self.store_outlier(read, reason)

        # get rid of node files.
        self.remove_node_files(node_id)

        # it is always sad to pop things
        self.nodes.pop(node_id)

        # and this recursion right here scares the shit out of me:
        if parent.node_id != 'root' and not parent.children:
            self.remove_node(parent.node_id)

    def store_outlier(self, unique_read_object, reason='unknown_reason'):
        if reason not in self.outlier_reasons:
            self.outlier_reasons.append(reason)
            self.outliers[reason] = set([])

        self.outliers[reason].add(unique_read_object)

    def store_outliers(self, unique_read_objects, reason='unknown_reason'):
        if reason not in self.outlier_reasons:
            self.outlier_reasons.append(reason)
            self.outliers[reason] = set([])

        self.outliers[reason].update(unique_read_objects)

    def remove_node_files(self, node_id):
        node = self.nodes[node_id]

        try:
            os.remove(node.alignment)
            os.remove(node.entropy_file)
            os.remove(node.unique_alignment)
        except:
            pass

    def merge_nodes(self, absorber_node_id, absorbed_node_id):
        # absorbed node gets merged into the absorber node
        absorber = self.get_node(absorber_node_id)
        absorbed = self.get_node(absorbed_node_id)

        absorber_parent = self.get_node(absorber.parent)
        absorbed_parent = self.get_node(absorbed.parent)

        # append read objects of absorbed to absorber:
        absorber.reads += absorbed.reads
        absorber.dirty = True

        # remove absorbed from the topology
        self.logger.info('topology: remove child "%s" from parent "%s"' % (absorbed.node_id, absorbed_parent.node_id))
        absorbed_parent.children.remove(absorbed.node_id)

        if absorber_parent.node_id != absorbed_parent.node_id and absorbed_parent.node_id != 'root':
            absorbed_parent.size -= absorbed.size

        # did the absorbed node's parent just become a final node?
        #  if that's the case we gotta remove this from the topology, because reads in this intermediate
        #  node was previously split between its child nodes. if all children were absorbed by other nodes
        # this node has no place in the topology anymore.
        if not absorbed_parent.children:
            self.logger.info('topology: remove empty parent of absorbed "%s"' % (absorbed_parent.node_id))
            self.remove_node(absorbed_parent.node_id)

        self.final_nodes.remove(absorbed.node_id)
        self.alive_nodes.remove(absorbed.node_id)

        #  if absorbed was a identified as a zombie:
        if absorbed.node_id in self.zombie_nodes:
            self.zombie_nodes.remove(absorbed.node_id)

        # kill the absorbed
        absorbed.killed = True
        self.remove_node_files(absorbed.node_id)

        #  refresh the absorbing node
        absorber.refresh()

    def get_parent_node(self, node_id):
        return self.nodes[self.nodes[node_id].parent]

    def update_final_nodes(self):
        self.alive_nodes = [n for n in sorted(self.nodes.keys()) if not self.nodes[n].killed]

        # get final nodes sorted by abundance        
        final_nodes_tpls = [(self.nodes[n].size, n) for n in self.alive_nodes if not self.nodes[n].children]
        final_nodes_tpls.sort(reverse=True)
        self.final_nodes = [n[1] for n in final_nodes_tpls]

    def store_node_representatives(self, node_ids, output_file_path, store_gaps=False):
        output = u.FastaOutput(output_file_path)
        for node_id in node_ids:
            output.write_id(node_id)
            if store_gaps:
                output.write_seq(self.nodes[node_id].representative_seq, split=False)
            else:
                output.write_seq(self.nodes[node_id].representative_seq.replace(b'-', b''), split=False)
        output.close()

    def store_final_nodes(self):
        for node_id in self.final_nodes:
            node = self.get_node(node_id)
            node.store()

    def recompute_nodes(self):
        for node_id in self.final_nodes:
            node = self.get_node(node_id)
            if node.dirty:
                node.refresh()

    def get_best_matching_node(self, sequence, distance_node_tuples):
        #  here we have a sequence, and a number of nodes that was previously decided to be OK
        #  candidates for this read. this function will do the trick to find really the best one
        #  among all, even if it contradicts with the distance (rep seq of one node can be more
        #  distant from the read than another one, but the distant one may be a better node to put
        # the outlier read because it may be more appropriate for biological reasons, so it is a
        # good idea to align them, see whether differences coincide with discriminant locations
        # that were originally used for decomposition, etc)

        # FIXME: now we return the smallest distance, but it should get smarter at some point:
        distance_node_tuples.sort()
        return distance_node_tuples[0][1]

    def relocate_outlier(self, outlier_read_object, target_node_id, original_removal_reason):
        '''
            add an outlier read to an existing node (target_node_id).
                                                                        '''
        node = self.nodes[target_node_id]

        #  update node alignment file with the outlier
        node.reads.append(outlier_read_object)
        node.dirty = True

        # remove outlier from outliers object
        self.outliers[original_removal_reason].remove(outlier_read_object)


class Node:
    # def __init__(self, node_id, output_directory):
    def __init__(self, node_id, output_directory, logger):  # agc
        self.logger = logger  # agc
        self.node_id = node_id
        self.pretty_id = None
        self.reads = []
        self.representative_seq = None
        self.killed = False
        self.dirty = False
        self.entropy = None
        self.entropy_tpls = None
        # normalized_m is a heuristic to get away from the problem of
        # having very low abundance nodes to be stuck in nodes of
        # extreme abundance. normalization will take the most abundant
        # unique sequence in the dataset for normalization. 
        self.normalized_m = None
        self.parent = None
        self.children = []
        self.discriminants = None
        self.max_entropy = None
        self.average_entropy = None
        self.size = 0
        self.level = None
        self.density = None
        self.freq_curve_img_path = None
        self.competing_unique_sequences_ratio = None
        self.file_path_prefix = os.path.join(output_directory, node_id)
        self.alignment_path = self.file_path_prefix + '.fa'
        self.unique_alignment_path = self.file_path_prefix + '.unique'

    def __str__(self):
        return self.pretty_id

    def set_representative(self):
        self.reads.sort(key=lambda x: x.frequency, reverse=True)
        self.representative_seq = self.reads[0].seq

    def do_entropy(self):
        self.entropy_tpls = []
        for position in range(0, len(self.representative_seq)):
            column = ''.join([read.seq.encode()[position] * read.frequency for read in self.reads])

            if len(set(column)) == 1:
                self.entropy_tpls.append((position, 0.0), )
            else:
                e = entropy(column)

                if e < 0.00001:
                    self.entropy_tpls.append((position, 0.0), )
                else:
                    self.entropy_tpls.append((position, e), )

        self.entropy = [t[1] for t in self.entropy_tpls]
        # self.logger.info(self.entropy)  # agc
        self.entropy_tpls = sorted(self.entropy_tpls, key=operator.itemgetter(1), reverse=True)
        self.max_entropy = max(self.entropy)
        self.average_entropy = numpy.mean([e for e in self.entropy if e > 0.05] or [0])

    def do_entropy_agc(self):  # agc
        self.logger.info("calculate entropy for node %s with size %d..." % (self.pretty_id, self.size))
        multiplex_reads_seq = numpy.empty([self.size, len(self.representative_seq)], dtype=numpy.ubyte)
        it = -1
        for i, read in enumerate(self.reads):
            # list_seq = list(read.seq)
            list_seq = bytearray(read.seq)
            for j in range(read.frequency):
                it += 1
                multiplex_reads_seq[it, :] = list_seq

        # 65, 71, 84, 67, 45  ==> AGTC-
        p = numpy.empty([5, len(self.representative_seq)], dtype=numpy.float64)
        test1 = multiplex_reads_seq == 65
        half_test = False
        if isinstance(test1, bool):
            # too large array try half
            half_test = True
            gc.collect()

        if half_test:
            test2_1 = multiplex_reads_seq[:self.size//2, :] == 65
            if isinstance(test2_1, bool):
                self.logger.info("Not enough memory to calculate entropy with matrix method")
                del multiplex_reads_seq
                self.do_entropy()
                return

            p[0, :] = numpy.count_nonzero(test2_1, axis=0) + \
                      numpy.count_nonzero(multiplex_reads_seq[self.size//2:, :] == 65, axis=0)
            del test2_1
            p[1, :] = numpy.count_nonzero(multiplex_reads_seq[:self.size//2, :] == 84, axis=0) +\
                      numpy.count_nonzero(multiplex_reads_seq[self.size//2:, :] == 84, axis=0)
            p[2, :] = numpy.count_nonzero(multiplex_reads_seq[:self.size//2, :] == 67, axis=0) + \
                      numpy.count_nonzero(multiplex_reads_seq[self.size // 2:, :] == 67, axis=0)
            p[3, :] = numpy.count_nonzero(multiplex_reads_seq[:self.size//2:, :] == 71, axis=0) + \
                      numpy.count_nonzero(multiplex_reads_seq[self.size // 2:, :] == 71, axis=0)
            p[4, :] = numpy.count_nonzero(multiplex_reads_seq[:self.size//2:, :] == 45, axis=0)+ \
                      numpy.count_nonzero(multiplex_reads_seq[self.size // 2:, :] == 45, axis=0)
        else:
            p[0, :] = numpy.count_nonzero(test1, axis=0)
            p[1, :] = numpy.count_nonzero(multiplex_reads_seq == 84, axis=0)
            p[2, :] = numpy.count_nonzero(multiplex_reads_seq == 67, axis=0)
            p[3, :] = numpy.count_nonzero(multiplex_reads_seq == 71, axis=0)
            p[4, :] = numpy.count_nonzero(multiplex_reads_seq == 45, axis=0)

        # self.logger.info("Counts A %s" % "".join(["%d " % p_ for p_ in p[0, :]]))
        # self.logger.info("Counts G %s" % "".join(["%d " % p_ for p_ in p[1, :]]))
        # self.logger.info("Counts T %s" % "".join(["%d " % p_ for p_ in p[2, :]]))
        # self.logger.info("Counts C %s" % "".join(["%d " % p_ for p_ in p[3, :]]))
        # self.logger.info("Counts - %s" % "".join(["%d " % p_ for p_ in p[4, :]]))

        p = p / self.size
        # self.logger.info("Size - %d" % self.size)
        # self.logger.info("Prob A %s" % "".join(["%g " % p_ for p_ in p[0, :]]))
        # self.logger.info("Prob G %s" % "".join(["%g " % p_ for p_ in p[1, :]]))
        # self.logger.info("Prob T %s" % "".join(["%g " % p_ for p_ in p[2, :]]))
        # self.logger.info("Prob C %s" % "".join(["%g " % p_ for p_ in p[3, :]]))
        # self.logger.info("Prob - %s" % "".join(["%g " % p_ for p_ in p[4, :]]))

        #  -(p_c*np.log(p_c) + p_a*np.log(p_a) + p_b*np.log(p_b) + p_t*np.log(p_t))
        p = numpy.multiply(p, numpy.ma.log2(p).filled(0))
        p[numpy.isnan(p) | numpy.isinf(p)] = 0

        # self.logger.info("ent A %s" % "".join(["%g " % p_ for p_ in p[0, :]]))
        # self.logger.info("ent G %s" % "".join(["%g " % p_ for p_ in p[1, :]]))
        # self.logger.info("ent T %s" % "".join(["%g " % p_ for p_ in p[2, :]]))
        # self.logger.info("ent C %s" % "".join(["%g " % p_ for p_ in p[3, :]]))
        # self.logger.info("ent - %s" % "".join(["%g " % p_ for p_ in p[4, :]]))

        p = -numpy.sum(p, axis=0)
        # Ignore too small numbers
        p = numpy.where(p < 0.00001, 0.0, p)

        self.entropy = p.tolist()
        self.entropy_tpls = [(i, v) for (i, v) in enumerate(p)]
        self.entropy_tpls = sorted(self.entropy_tpls, key=operator.itemgetter(1), reverse=True)
        self.max_entropy = max(self.entropy)
        self.average_entropy = numpy.mean([e for e in self.entropy if e > 0.05] or [0])
        # log the entropy 
        # self.logger.info("Entropy %s" % "".join(["%g " % en for en in self.entropy]))
        self.logger.info("Max entropy %g" % self.max_entropy)

    def set_normalized_m(self, default_min_entropy, most_abundant_unique_sequence_in_the_dataset):
        """When called, this function sets the normalized_m value for the node.
           Normalized m heuristic can be used for decomposition of a given node
           instead of the static m that is provided by the user"""

        y = numpy.sqrt(self.size) * 1.0 / numpy.sqrt(most_abundant_unique_sequence_in_the_dataset * 2)
        self.normalized_m = default_min_entropy - (default_min_entropy * y)
        if self.normalized_m < 0.05:
            self.normalized_m = 0.05

    def do_competing_unique_sequences_ratio_and_density(self):
        if len(self.reads) == 1:
            self.competing_unique_sequences_ratio = 0
        else:
            self.competing_unique_sequences_ratio = self.reads[1].frequency * 1.0 / self.reads[0].frequency

        self.density = self.reads[0].frequency * 1.0 / self.size

    def refresh(self):
        self.set_representative()
        self.size = sum([read.frequency for read in self.reads])

        # self.do_entropy()  # agc
        self.do_entropy_agc()

        self.do_competing_unique_sequences_ratio_and_density()
        self.dirty = False

    def store(self):
        alignment = open(self.alignment_path, 'w')
        for read in self.reads:
            for read_id in read.ids:
                alignment.write('>%s\n%s\n' % (read_id, read.seq.decode()))
        alignment.close()

        unique_alignment = open(self.unique_alignment_path, 'w')
        for read in self.reads:
            unique_alignment.write('>%s|frequency:%d\n%s\n' % (read.ids[0], read.frequency, read.seq.decode()))
        unique_alignment.close()
