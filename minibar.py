#!/usr/bin/env python

from __future__ import division, print_function
import sys, errno, timeit

try:
    import edlib
except:
    print("\n    edlib module not found, to install use:\n\n        pip install edlib\n", file=sys.stderr)
    exit(3)


def read_barcode_file(primer_filename, ops):
    global fwd_indexes, rev_indexes
    global len_first_index
    global fwd_primer, len_fwd_primer, rev_primer, len_rev_primer
    global indx2sample_map_fr, indx2sample_map_rf
    global fwdonly_sampleID # for method 0 which only requires a single forward barcode

    # read the primer file contents into a list
    def read_primer_list(primer_filename):
        def header_line(ln, ops):
            if ops.first_line_header == -1:  # -1 means auto detect
                lwln = ln.lower()
                return lwln.find('sample') >= 0 or lwln.find('id') >= 0
            else:
                return ops.first_line_header != 0

        primers = []
        headers = []
        with open(primer_filename) as primerFile:
            ln1 = primerFile.readline().rstrip()
            ln1_list = ln1.split("\t")
            if header_line(ln1, ops):
                headers = ln1_list
            else:
                primers.append(ln1_list)

            if len(ln1_list) < 5:
                print("Need at least 5 tab delimited columns in the barcode_file.\n"
                      "Here is the first line of '{}':\n{}\n".format(primer_filename, ln1), file=sys.stderr)
                sys.exit(2)

            for line in primerFile:
                primers.append(line.rstrip().split("\t"))
        return primers

    def set_col_ix(primer_list, ops):
        if ops.col_set_on_cmdline or len(primer_list) < 1:
            return

        # if user didn't set explicitly, set col pos based on number of fields
        ops.sampleix = 0
        num_cols = len(primer_list[0])
        if num_cols == 5:
            ops.fwix = 1; ops.fwprm = 2; ops.rvix = 3; ops.rvprm = 4
        elif num_cols == 6:  # presume 1 is a name
            ops.fwix = 2; ops.fwprm = 3; ops.rvix = 4; ops.rvprm = 5
        elif num_cols in [7,8,9]:  # presume 1 and 4 are names
            ops.fwix = 2; ops.fwprm = 3; ops.rvix = 5; ops.rvprm = 6
        else:  # leave as the 11 column defaults: SAMPLE = 0; FWIX = 2; FWPRM = 4; RVIX = 7; RVPRM = 9
            return #ops.fwix = 2; ops.fwprm = 4; ops.rvix = 7; ops.rvprm = 9

    # make 2 dict maps, keys are a combination of the two indexes that
    # uniquely define the Sample. we use 2 maps so we don't care if the
    # first index found is the forward or reverse index.
    def make_sample_id_maps(primer_list, sampleid_col, fwix_col, rvix_col):
        fwixmap = {}; rvixmap = {}; fwdonly_ids = {}
        for p in primer_list:
            id = p[sampleid_col].strip()
            if id == "": continue
            fwix = p[fwix_col].upper(); rvix = p[rvix_col].upper()
            fwixmap[fwix + '|' + rvix] = id
            rvixmap[rvix + '|' + fwix] = id
            if not fwix in fwdonly_ids: # for method 0 remember the first sample associated with the forward index
                fwdonly_ids[fwix] = id
        return fwixmap, rvixmap, fwdonly_ids

    # (fwd and rv) index will occur multiple times in the primer_list
    # this will make a dictionary entry for each unique index found at ix_col.
    # each entry lists the indices in primer_list where this index occurs
    # also UPPERCASE the index since sequence will be in uppercase
    def make_index_list(primer_list, ix_col):
        index_map = {}
        for i, p in enumerate(primer_list):
            index = p[ix_col].upper()
            if index in index_map:
                index_map[index].append(i)
            else:  # first time seeing this index
                index_map[index] = [i]
        return index_map
    
    def check_method_0_dups():
        if not ops.search_method == 0:
            return
        dup = 0
        for fwix in fwdonly_sampleID:
            if len(fwd_indexes[fwix]) > 1:
                if dup==0: print("", file=sys.stderr)
                dup += 1
                print("Index '{}' occurs with multiple samples. Method 0 output is for sample '{}'.".format(fwix, fwdonly_sampleID[fwix]), file=sys.stderr)
        if dup: print("", file=sys.stderr)

    # 02Dec2018 JBH v0.21 use of same index in fwd and reverse won't work for Method 2
    def check_dup_fwd_rev_index(fwd_indexes, rev_indexes):
        dup_index = []
        for f in fwd_indexes:
            for r in rev_indexes:
                if f.upper() == r.upper():
                    dup_index.append(f)
        if ops.search_method != 1 and len(dup_index) > 0:
            errstr = "Using same index in forward and reverse lists: " + " ".join(dup_index)
            report_duplicates = error if ops.search_method == 2 else warning
            if ops.search_method == 2:
                errstr = "Not valid for Method 2: " + errstr
            report_duplicates(errstr)

    # make sure ids and index/pairs not blank and not duplicates; and number of cols is the same for each line
    def check_empty_and_dups(primer_list, sampleid_col, fwix_col, rvix_col):
        sample_map = {}; index_pair_map = {}; col_count_map = {}; col_count_first_line = {}
        missing_samples = []; dup_samples = []; dup_indexpairs = []
        for i, p in enumerate(primer_list):
            id = p[sampleid_col].strip()
            if id == "":
                missing_samples.append(i)
            elif id in sample_map:
                dup_samples.append([id, i])
            else:
                sample_map[id] = i

            if fwix_col >= len(p) or rvix_col >= len(p):  # not enough tabbed fields in line
                error("Not enough tabbed fields in entry {}: {}".format(i+1,p))

            fwix = p[fwix_col].upper(); rvix = p[rvix_col].upper()
            ix_pair = fwix + ' | ' + rvix
            if ix_pair in index_pair_map:
                dup_indexpairs.append([ix_pair, i])
            else:
                index_pair_map[ix_pair] = id

            col_count = len(p)
            if not col_count in col_count_map:
                col_count_map[col_count] = 0
                col_count_first_line[col_count] = "\t".join(p)
            col_count_map[col_count] += 1

        #  report any errors
        errstr = ""
        if len(missing_samples) > 0:
            errstr = "No sample names on items: {}\n    ".format(", ".join(map(str, missing_samples)))
        for dup in dup_samples:
            errstr += "{} sample name duplicated on item {}\n    ".format(dup[0], dup[1])
        for dup in dup_indexpairs:
            errstr += "{} index pairs duplicated on item {}\n    ".format(dup[0], dup[1])
        if len(col_count_map) > 1:
            expected = len(primer_list[0])  # number of columns we expected in each line
            errstr += "Some lines in the barcode file have other than {} columns:\n".format(expected)
            example_str = "    {} line(s) with {} columns, e.g: {}\n"
            for c in col_count_map:
                if c != expected:
                    errstr += example_str.format(col_count_map[c], c, col_count_first_line[c])

        if errstr != "":
            report_duplicates = error if ops.error_on_duplicate_samples else warning
            report_duplicates(errstr)

    # read the tsv data from the primer and index file
    primers = read_primer_list(primer_filename)

    # define column indices for the Index and the Primer both Forward and Backward
    set_col_ix(primers, ops) # set based on number of columns if user didn't specifically tell us what to use
    FWIX = ops.fwix; RVIX = ops.rvix; FWPRM = ops.fwprm; RVPRM = ops.rvprm; SAMPLE = ops.sampleix

    # do some basic error checking
    check_empty_and_dups(primers, SAMPLE, FWIX, RVIX)

    # gather unique set of fwd and rev indexes
    fwd_indexes = make_index_list(primers, FWIX)
    rev_indexes = make_index_list(primers, RVIX)
    first_index = primers[0][FWIX]
    len_first_index = len(first_index)

    # 02Dec2018 JBH check to see if same index/barcode used in fwd_indexes and rev_indexes, issue warning
    check_dup_fwd_rev_index(fwd_indexes, rev_indexes)

    # set fwd and rev primer vars
    fwd_primer = primers[0][FWPRM]
    len_fwd_primer = len(fwd_primer)
    rev_primer = primers[0][RVPRM]
    len_rev_primer = len(rev_primer)

    # map index pairs to sample_id. create 2 dicts one with FWIX|RVIX key other with RVIX|FWIX key
    indx2sample_map_fr, indx2sample_map_rf, fwdonly_sampleID = make_sample_id_maps(primers, SAMPLE, FWIX, RVIX)
    check_method_0_dups() # 09Apr2020 flag as ambiguous dup fwd barcodes/indexes if using Method
    return primers


def display_barcode_file_inf(ops):
    if len(ops.barcode_file_info) < 1 or not ops.barcode_file_info[0] in 'frpcba':
        error('"{}" is an invalid argument for -info'.format(ops.barcode_file_info))

    def display_col_info(primer_list, ops):
        first = primer_list[0]
        cols = "Sample name     in col {}:\t{}\n".format(1+ops.sampleix, first[ops.sampleix])
        cols += "Forward Barcode in col {}:\t{}\n".format(1+ops.fwix, first[ops.fwix])
        cols += "Forward Primer  in col {}:\t{}\n".format(1+ops.fwprm, first[ops.fwprm])
        cols += "Reverse barcode in col {}:\t{}\n".format(1+ops.rvix, first[ops.rvix])
        cols += "Reverse Primer  in col {}:\t{}\n".format(1+ops.rvprm, first[ops.rvprm])
        print(cols)

    def display_index_editdistances(indexes):
        shortest = len_first_index+1; ixlen = len_first_index
        ln = 'Indexes '
        for i in indexes: # header with each index in a column
            ln += '\t' + i
        print(ln)

        for i in indexes:
            ln = i
            for j in indexes:
                result = edlib.align(i,j, task='distance')
                dist = result['editDistance']
                if dist > 0 and dist < shortest:
                    shortest = dist
                    ixlen = len(j)
                ln += '\t' + str(dist)
            print(ln)
        print('\nClosest has edit distance of {} with index lengths of {}, {:.2%} alike'.format(shortest, ixlen, 1-(shortest/ixlen)))

    primer_list = read_barcode_file(ops.primerfile, ops)

    info = ops.barcode_file_info[0]
    if info == 'p' or info == 'a':
        result = edlib.align(fwd_primer, rev_primer, task='path', additionalEqualities=IUPAC_maps)
        dist = result['editDistance']
        alike = 1-(dist/len_fwd_primer)
        msg = 'Primers {} {} edit distance {} with lengths of {}, {:.2%} alike'
        print(msg.format(fwd_primer, rev_primer, dist, len_rev_primer, alike))
    elif info == 'c':  # show cols values for first line 14Jun2018
        display_col_info(primer_list, ops)
    if info in 'frba':
        if info != 'b' and info != 'a':
            indexes = fwd_indexes if info != 'r' else rev_indexes
        else: # we want 'b'oth indexes combined to see what we have
            indexes = fwd_indexes.copy()
            indexes.update(rev_indexes)
        display_index_editdistances(indexes)


def open_sequence_file(fqname):  # won't allow stdin with this method
    global fa_headerln
    if fqname[-3:] == ".gz":  # open as a gzipped file
        import gzip
        fh = gzip.open(fqname, 'rt')
    else:  # use regular open
        fh = open(fqname)

    firstln = fh.readline().rstrip()
    if firstln[0] == '@':  # it's a fastq file
        fh.seek(0)  # rewind file to beginning
        isfastq = True
    else:
        isfastq = False
        if firstln[0] == '>':  # it's a fasta file either
            fa_headerln = firstln
        else:
            print("{} is not a Fastq or a Fasta file.".format(fqname))
            return None, False

    global load_seqrecs
    load_seqrecs = load_fq_seqrecs if isfastq else load_fa_seqrecs

    return fh, isfastq


def load_fq_seqrecs(fh_seq, maxseqs):
    global rec_seqs, rec_hdrs, rec_quals  # appends to these (caller is responsible for emptying the lists)

    RECSIZE = 4
    rec_ln = 0
    for ln in fh_seq:
        ln = ln.rstrip()
        rec_ln += 1
        if rec_ln == 1:  # header line
            rec_hdrs.append(ln)
        elif rec_ln == 2:
            rec_seqs.append(ln)
        elif rec_ln == 4:
            rec_quals.append(ln)
        if rec_ln == RECSIZE:
            rec_ln = 0
            if len(rec_seqs) >= maxseqs:
                break

    return len(rec_seqs)


def load_fa_seqrecs(fh_seq, maxseqs):
    global rec_seqs, rec_hdrs, rec_quals  # appends to these (caller is responsible for emptying the lists)
    global fa_headerln  # header line for next fasta record

    faseq = ''
    for ln in fh_seq:
        ln = ln.rstrip()
        if ln[0] != '>':
            faseq += ln
        else: # got to header of next fasta record
            rec_hdrs.append(fa_headerln)
            rec_seqs.append(faseq)
            rec_quals.append('') # placeholder for fasta, not used

            fa_headerln = ln
            faseq = ''
            if len(rec_seqs) >= maxseqs:
                break

    if faseq != '': # presume it is because we reached end of file, not maxseqs
        rec_hdrs.append(fa_headerln)
        rec_seqs.append(faseq)
        rec_quals.append('') # placeholder for fasta, not used

    return len(rec_seqs)


def rev_comp_py2(nt_str):
    import string
    trans_tbl = string.maketrans('AaTtGgCcNn', 'TtAaCcGgNn')
    return nt_str.translate(trans_tbl)[::-1]


def rev_comp_py3(nt_str):
    trans_tbl = str.maketrans('AaTtGgCcNn', 'TtAaCcGgNn')
    return nt_str.translate(trans_tbl)[::-1]


# v 0.20 20Nov2018 add 3rd arg to determine which set (reverse or forward) to check first
# this will support using same index pairs but flipped from forward to reverse to id a sample
def sample_id_from_indexes(index1, index2, check_rev_first=None):
    key = index1 + '|' + index2
    map_one = indx2sample_map_fr
    map_two = indx2sample_map_rf
    if check_rev_first:
        map_one, map_two = map_two, map_one

    if key in map_one:
        return map_one[key]
    elif key in map_two:
        return map_two[key]
    return ''

#    if key in indx2sample_map_fr:
#        return indx2sample_map_fr[key]
#    elif key in indx2sample_map_rf:
#        return indx2sample_map_rf[key]
#    return ''


def ids_from_index_matches(ind_match_1, ind_match_2):
    id_list = []; id_dict = {}
    for i1 in ind_match_1:
        for i2 in ind_match_2:
            id = sample_id_from_indexes(i1[0], i2[0])
            if id != '' and not id in id_dict:
                id_list.append(id)
                id_dict[id] = True

    return id_list


# 06Dec2022 add N and other IUPAC codes that might be in a primer to YR, 2 base ones WSMK
IUPAC_maps = [("Y", "C"), ("Y", "T"), ("R", "A"), ("R", "G"),
              ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T"),
              ("W", "A"), ("W", "T"), ("M", "A"), ("M", "C"),
              ("S", "C"), ("S", "G"), ("K", "G"), ("K", "T")]

def primer_positions(seq, primer):
    rs = edlib.align(primer, seq, 'HW', 'locations', max_dist_primer, additionalEqualities=IUPAC_maps)
    return rs['editDistance'], rs['locations']


# this uses edlib.align()
max_dist_index = 4
max_search_area = 80
def search_seq_for_indexes(seq, indexes):
    rslts = []
    seq_prefix = seq[:max_search_area]
    for query in indexes:
        rs = edlib.align(query, seq_prefix, 'HW', 'locations', max_dist_index)
        dist = rs['editDistance']
        if dist > -1:
            rslts.append([query, dist, rs['locations']])
    return rslts


# index_matches are zero, one or more matches of an index query
#   (if more than 1, they will have have the same dist score).
# primer_matches first element is the dist score, then a list of
#   matches with that same score, each entry being a tuple of beg/end pos,eg (47, 73)
#
# best case is an index match with a good primer match right after it.
# return is a tuple of best match and primer best match, if any.
# good hit ex: (1, 'CGCTCTGCCAAAGAT', (33, 46), 3, (47, 73))
# where we have index dist score, index, pos, primer dist score, primer pos
# -1, '', () for index no-match; -1, () for primer no match
def choose_best_index(index_matches, primer_matches):
    CLOSE_ENOUGH = 5
    best_score = 999; best_im = None; best_loc = None; primer_loc = None
    no_result = (-1, '', (), -1, ())

    if len(index_matches) == 0 or len(primer_matches) == 0:
        return no_result

    primer_score = primer_matches[0]
    if primer_score > -1:
        # loop over primers first, since there is usually 1 or none rather than more, unlike index close matches
        for pm in primer_matches[1]:  # each pm is a tuple of the primer beg, end position
            pm_beg = pm[0]; closest = 999
            for im in index_matches:
                score = im[1]
                for loc in im[2]:
                    close = abs(abs(pm_beg - loc[1])-1)
                    if close < CLOSE_ENOUGH:  # good primer match right after good index match
                        if score < best_score or (score==best_score and close < closest):  # better score or same score and closer
                            best_im = im
                            best_loc = loc
                            primer_loc = pm
                            closest = close
                            best_score = score

    return (best_im[1], best_im[0], best_loc, primer_score, primer_loc) if best_im else no_result


def find_best_index(seq, fwd=True):
    global ind_matches, prm_matches
    indexes = fwd_indexes if fwd else rev_indexes
    primer = fwd_primer if fwd else rev_primer

    ind_matches = search_seq_for_indexes(seq, indexes)
    prm_matches = primer_positions(seq, primer)
    return choose_best_index(ind_matches, prm_matches)


def search_for_best_index(seq, skip_fwdchk=False):
    global best_index

    strand = '?'; best = [-1]
    if not skip_fwdchk:
        best = find_best_index(seq, True)  # true means use fwd indexes

    if best[0] == -1:  # try rev indexes and primer
        best = find_best_index(seq, False)  # False means use rev indexes
        if best[0] != -1:
            best_index = best[1]
            strand = '-'
    else:  # found it with fwd indexes
        best_index = best[1]
        strand = '+'

    return (strand,) + best


# for each sequence in the rec_seqs this will output results for the sample ID based
# on finding a good fwd_index/fwd_primer or rev_index/rev_primer pair at the beginning
# of the sequence. the complementary pair at the other end is best find, but if that
# is not there then we use complementary matched indexes to determine sample IDs
def search_sequence_list(seq_number, rec_seqs, rec_hdrs, rec_quals, ops, fh_map):
    global H, HH, Hh, hh, hx, samps, mult_ids

    USE_FWD_INDEXES = True; USE_REV_INDEXES = False
    HH_method = 1; hh_method = 2; Hh_method = 3; single_h_method = 0  # add 0 method 31Mar2020
    unknown_id = 'unk'
    outtype = ops.output
    search_method = ops.search_method
    assert(search_method in [HH_method, Hh_method, hh_method, single_h_method])
    show_color = ops.show_color

    seq_ix = -1
    for fseq in rec_seqs:
        seq_ix += 1; seq_number += 1

        ID_matches = 0; all_ids = ''; strength = '  '
        ix2 = ''
        ix1_loc = prm1_loc = ix2_loc = prm2_loc = rslt = None

        strand = '?'  # easy way to shoehorn in new method without deeper indentation
        if not search_method in [hh_method, single_h_method]:   # replaced search_method != hh_method
            rslt = search_for_best_index(fseq)  # looking for a good match of an index close to a primer
            strand = rslt[0]

        if strand != '?':  # strong hit for beginning index & primer, check ending index
            H += 1
            ix1 = best_index
            ixed1 = rslt[1]; prmed1=rslt[4]; ed1='('+str(ixed1)+','+str(prmed1)+')'
            ix1_loc = rslt[3]; prm1_loc = rslt[5]
            rc_fseq = rev_comp(fseq)
            if strand == '+':  # found fwd indexes at beginning of seq, search revcomp of seq for rev indexes
                strand2 = '-'
                best = find_best_index(rc_fseq, USE_REV_INDEXES)
            else:  # found rev indexes at beginning of seq, search revcomp of seq for fwd indexes
                assert strand == '-'
                strand2 = '+'
                best = find_best_index(rc_fseq, USE_FWD_INDEXES)

            if best[0] != -1:  # strong hit for ending index too
                rcrslt = (strand2,) + best
                strength = 'HH'
                HH += 1
                ix2 = best[1]
                ix2_loc = rcrslt[3]; prm2_loc = rcrslt[5]
                ed2 = '('+str(best[0])+','+str(best[3])+')'
            else:  # no primer hit, see if there are any index matches
                rcrslt = (strand2,) + tuple(ind_matches)
                if len(ind_matches) == 0:
                    strength = 'Hx'
                    ed2 = '(-1,-1)'
                else:
                    Hh += 1
                    # loop through matches and choose those that create a sample id
                    best_score = 999; best_im = None; ID_matches = 0
                    for im in ind_matches:
                        score = im[1]
                        if score <= best_score:  # if none yet, or none that matches a sample, choose this im
                            best_score = score
                        if not best_im:
                            best_im = im
                        id = sample_id_from_indexes(ix1, im[0], strand == '-')
                        if id != '':
                            ID_matches += 1
                            all_ids += id + ' '
                            if ID_matches == 1 or score < best_score:
                                # all_ids = id
                                best_im = im  # take first or higher scoring match as best
                    ix2 = best_im[0]  # sequence
                    ix2_loc = best_im[2][0]
                    strength = 'Hh'
                    ed2 = '({},-1)'.format(best_im[1])  # score
            sample_id = sample_id_from_indexes(ix1, ix2, strand == '-') if ID_matches <= 1 else all_ids.strip()
            if sample_id != '': samps += 1
            else: sample_id = unknown_id
            if ID_matches > 1: mult_ids += 1
        elif search_method in [hh_method, Hh_method, single_h_method]: # replaced search_method == hh_method or search_method == Hh_method
            strand2 = '?'
            ed1 = '(-1)'; ed2 = '(-1)'
            rcrslt = (-1, '', (), -1, ())

            rc_fseq = rev_comp(fseq)

            beg_fwd = search_seq_for_indexes(fseq, fwd_indexes)
            end_rev = search_seq_for_indexes(rc_fseq, rev_indexes)
            fr_ids = ids_from_index_matches(beg_fwd, end_rev)
            if len(fr_ids) > 0:
                strand  = '+'
                strand2 = '-'
                rslt = (strand,) + tuple(beg_fwd)
                rcrslt = (strand2,) + tuple(end_rev)
                strength = 'hh'
                ed1 = '({})'.format(beg_fwd[0][1])
                ed2 = '({})'.format(end_rev[0][1])

            beg_rev = search_seq_for_indexes(fseq, rev_indexes)
            end_fwd = search_seq_for_indexes(rc_fseq, fwd_indexes)
            rf_ids = ids_from_index_matches(beg_rev, end_fwd)
            if len(rf_ids) > 0 and strand == '?':
                strand = '-'
                strand2 = '+'
                rslt = (strand,) + tuple(beg_rev)
                rcrslt = (strand2,) + tuple(end_fwd)
                strength = 'hh'
                ed1 = '({})'.format(beg_rev[0][1])
                ed2 = '({})'.format(end_fwd[0][1])

            all_id_list = fr_ids + rf_ids

            if len(all_id_list) > 0:
                samps += 1
                hh += 1
                sample_id = " ".join(all_id_list)
                if len(all_id_list) > 1:
                    mult_ids += 1
                    ID_matches = len(all_id_list) # 19Oct2022 was causing mult ID named files in this scenario
                ix1_loc = rslt[1][2][0]
                ix2_loc = rcrslt[1][2][0]
            elif search_method == single_h_method and beg_fwd: # hx or xh  # remove beg_rev for now
                # 26Mar2020 implement any port in a storm mode, which lets us grab onto any found index
                # min = a if a < b else b
                samps += 1
                beg_matched = beg_fwd if beg_fwd else beg_rev
                ixseq = beg_matched[0][0]
                sample_id = fwdonly_sampleID[ixseq] if ixseq in fwdonly_sampleID else "unk"
                hx += 1
                if beg_fwd:
                    strand = '+'; strand2 = '-'
                    strength = "hx"
                    ed1 = '({})'.format(beg_matched[0][1])
                else:
                    strand = '-'; strand2 = '+'
                    strength = "xh"
                    ed2 = '({})'.format(beg_matched[0][1])
                rslt = (strand,) + tuple(beg_matched)
                ix1_loc = rslt[1][2][0]
            else:
                strand = 'x'; strand2 = 'x'  # report total misses as x(-1), x(-1) not ?(-1), ?(-1)
                sample_id = unknown_id
        else:
            ed1 = '(-1,-1)'; ed2 = '(-1,-1)'
            if rslt and int(rslt[1]) > -1:
                ed1 = '({},{})'.format(rslt[1],rslt[2])
            strand2 = 'x'; strand = 'x'  # report total misses as x(-1,-1), x(-1,-1) not ?(-1,-1), ?(-1,-1)
            sample_id = unknown_id
            rcrslt = (-1, '', (), -1, ())

        if outtype == 3:  # 3 is the diagnostic output mode
            if strand != '?':
                print(seq_number, sample_id, '\t', strength, strand, rslt[1:], strand2, rcrslt[1:])
        else:
            if outtype == 2:  # use upper/lower case to distinguish index in sequence, and potentially color
                fseq = make_display_seq(fseq, ix1_loc, prm1_loc, ix2_loc, prm2_loc, show_color)
            elif outtype == 4:  # Trim primers from sequence and if fastq also trim qual line
                beg_seq, end_seq = get_trim_locs(fseq, ix1_loc, prm1_loc, ix2_loc, prm2_loc)
                fseq = fseq[beg_seq: end_seq]
                if ops.isfastq:
                    rec_quals[seq_ix] = rec_quals[seq_ix][beg_seq: end_seq]

            output_seq(sample_id, ID_matches, fh_map, ops,
                       fseq, rec_hdrs[seq_ix], rec_quals[seq_ix],
                       strength[0]+strand+ed1+','+strength[1]+strand2+ed2)

    return seq_number


def output_seq(sample_id, ID_matches, fh_map, ops, seq, hdr, quals, match_info):

    fh = sys.stdout
    if ops.output_to_files:
        sample_name = sample_id if ID_matches <= 1 else "Multiple_Matches"  # see where we are writing the output
        fh = get_sample_fh(sample_name, fh_map, ops.output_file_prefix, ops.isfastq)

    fh.writelines([hdr,' ', match_info, ' ', sample_id, '\n'])
    fh.writelines([seq, '\n'])
    if ops.isfastq:
        fh.writelines(['+\n', quals, '\n'])


def get_trim_locs(seq, ix1_loc, prm1_loc, ix2_loc, prm2_loc):
    ln = len(seq)
    beg_seq = 0
    end_seq = ln

    if prm1_loc:  # set begin to char after end of left primer
        beg_seq = prm1_loc[1]
    elif ix1_loc:  # set begin to char after end of index + length of primer
        beg_seq = ix1_loc[0] + len_first_index + len_fwd_primer

    if prm2_loc:
        end_seq = ln-1 - prm2_loc[1]
    elif ix2_loc:
        end_seq = ln-1 - ix2_loc[1] - len_rev_primer

    return beg_seq, end_seq


blue = "\033[0;34m"; green = "\033[0;32m";  red = "\033[0;31m"; NC = "\033[0m"  # No Color (reset)
# flip case to highlight index and primer features in sequence, and apply color to these areas if requested
def make_display_seq(seq, ix1_loc, prm1_loc, ix2_loc, prm2_loc, show_colors = False, invert2s = True):
    if not ix1_loc:
        return seq

    clr_ix1 = ''; clr_prm1 = ''; clr_ix2 = ''; clr_prm2 = ''; reset = NC if show_colors else ''

    # set values for 5' index and primer
    ix1_beg = ix1_loc[0]; ix1_end = ix1_beg + len_first_index

    if show_colors:
        clr_ix1 = blue

    if prm1_loc:
        prm1_beg = prm1_loc[0]
        if prm1_beg < ix1_end and prm1_beg >= ix1_beg:  # 18Jan2023 BUG FIX:
            ix1_end = prm1_beg  # if primer match overlaps index, reduce index end so it does not
        prm1_end = prm1_loc[1]
        if show_colors: clr_prm1 = green
    else:
        prm1_beg = ix1_end
        prm1_end = prm1_beg + len_fwd_primer
        if show_colors: clr_prm1 = red

    # set values for 3' primer and index, if no ix2_loc then no prm2_loc
    if ix2_loc:
        ix2_beg = ix2_loc[0]; ix2_end = ix2_loc[1]
        if show_colors: clr_ix2 = blue

        if prm2_loc:
            prm2_beg = prm2_loc[0]; prm2_end = prm2_loc[1]
            if show_colors: clr_prm2 = green
    else:
        ix2_beg = ix2_end = prm2_loc = None

    if invert2s:
        ln = len(seq)
        if ix2_loc:
            ib = ix2_beg
            ix2_beg = ln-1 - ix2_end
            ix2_end = ln - ib
            if prm2_loc:
                pb = prm2_loc[0]
                prm2_beg = ln-1 - prm2_loc[1]
                prm2_end = ln - pb
            else:  # set primer2 right before ix2_loc
                prm2_beg = ix2_beg - len_rev_primer
                prm2_end = ix2_beg
                if show_colors: clr_prm2 = red

    if not ix2_loc:  # no 3' index or primer
        display_seq = seq[0:ix1_beg].lower() + \
                      clr_ix1 + seq[ix1_beg:ix1_end].upper() + \
                      clr_prm1 + seq[prm1_beg:prm1_end].lower() + reset + seq[prm1_end:].upper()
    else:
        display_seq = seq[0:ix1_beg].lower() + \
                      clr_ix1 + seq[ix1_beg:ix1_end].upper() + \
                      clr_prm1 + seq[prm1_beg:prm1_end].lower() + reset + seq[prm1_end:prm2_beg].upper() + \
                      clr_prm2 + seq[prm2_beg:prm2_end].lower() + reset + \
                      clr_ix2 + seq[ix2_beg:ix2_end].upper() + reset + seq[ix2_end:].lower()
    return display_seq


def get_sample_fh(sampID, fh_map, prefix, isfastq):

    def make_sample_filename(sampID, prefix, isfastq):
        sample = "".join([x if x.isalnum() or x in "._-$#" else "_" for x in sampID])
        name = prefix + sample
        name += '.fastq' if isfastq else '.fasta'
        return name

    if sampID in fh_map:
        return fh_map[sampID]

    sample_filename = make_sample_filename(sampID, prefix, isfastq)
    try:
        sample_fh = open(sample_filename, 'w')
        fh_map[sampID] = sample_fh
        return sample_fh
    except:
        return sys.stdout


def minibar(ops):
    global H, HH, Hh, hh, hx, samps, mult_ids
    global rec_seqs, rec_hdrs, rec_quals
    global max_dist_index, max_dist_primer, max_search_area

    read_barcode_file(ops.primerfile, ops)

    fh, isfastq = open_sequence_file(ops.sequencefile)
    if not fh:
        print("Problem opening '{}'".format(ops.sequencefile), file=sys.stderr)
        return

    # holds file handles for sample files
    fh_map = {}
    ops.isfastq = isfastq

    max_search_area = ops.search_len
    ixed = ops.index_edit_distance
    max_dist_index  = ixed if (ixed >= 1) else len_first_index - int(len_first_index * ixed)
    max_dist_primer = int(len_fwd_primer * 0.3 + 3) if ops.primer_edit_distance==-1 else ops.primer_edit_distance

    info_msg = '{} {} : Index edit dist {}, Primer edit dist {}, Search Len {}, Search Method {}, Output Type {}'
    print(info_msg.format(ops.primerfile, ops.sequencefile, max_dist_index, max_dist_primer,
                                     max_search_area, ops.search_method, ops.output_letter), file=sys.stderr)

    start_time = timeit.default_timer()

    sequence_block_size = 10000

    last_seq_to_output = ops.num_seqs
    all_seqs = last_seq_to_output < 0
    seq_num = 0; throwaway = 0
    if ops.start_seq > 1:
        throwaway = ops.start_seq - 1
        last_seq_to_output += throwaway
        rec_seqs = []; rec_hdrs = []; rec_quals = []
        seq_num = load_seqrecs(fh, throwaway)

    H = 0; HH = 0; Hh = 0; hh = 0; hx = 0; samps = 0; mult_ids = 0

    while all_seqs or seq_num < last_seq_to_output:
        to_read = sequence_block_size if all_seqs else min(sequence_block_size, last_seq_to_output-seq_num)

        rec_seqs = []; rec_hdrs = []; rec_quals = []
        if load_seqrecs(fh, to_read) < 1:
            break

        seq_num = search_sequence_list(seq_num, rec_seqs, rec_hdrs, rec_quals, ops, fh_map)
        print(seq_num, end='\r', file=sys.stderr); sys.stderr.flush()

    elapsed = timeit.default_timer() - start_time

    if ops.search_method == 1:
        hits = 'H {} HH {} Hh {}'.format(H, HH, Hh)
    elif ops.search_method == 2:
        hits = 'hh {}'.format(hh)
    elif ops.search_method == 0:  # 09Apr2020 assigns based on finding forward index alone if nothing works
        hits = 'H {} HH {} Hh {} hh {} hx {}'.format(H, HH, Hh, hh, hx)
    else:  # method 3
        hits = 'H {} HH {} Hh {} hh {}'.format(H, HH, Hh, hh)

    result_str = '\r{} seqs: {} IDs {} Mult_IDs {} ({:0.4f}s)'.format(seq_num-throwaway, hits, samps, mult_ids, elapsed)
    print(result_str, file=sys.stderr)

    for fh in fh_map:  # close the sample files if we had any open
        try:
            fh.close()
        except:
            pass


def version():
    # 0.25 18Jan2023 fix -C -CC problem when fwd index match shorter than actual index and overlaps with primer
    # 0.24 06Dec2022 add N and other IUPAC codes that might be in a primer to YR, 2 base ones WSMK
    # 0.23 19Oct2022 fix problem with multiple sample named files
    # 0.22 31Mar2020 add 0 method to identify sample from even a forward barcode (any port in a storm)
    # 0.21 02Dec2018 -info both added, for same index in fwd and rev lists: err Method 2, warn Method 3, silent Method 1
    # 0.20 20Nov2018 sample_id_from_indexes() 3rd arg to swap indexes btw fwd/rev so same index in rev fwd works
    # 0.19 14Jun2018 -info cols
    # 0.18 12Jun2018 remove single, tighten file read # 0.17 10Jun2018 -T -w added
    return "minibar.py version 0.25"


def error(errmsg, exit_code=3):
    print('\n    ' + errmsg + '\n', file=sys.stderr)
    sys.exit(exit_code)


def warning(warnmsg):
    print('\n    ' + warnmsg + '\n', file=sys.stderr)


def usage(show_all_descrips=False):
    secondary_output_formats = \
        "-D diagnostic output, instead of sequence displays edit distances of index and primer matches\n"

    secondary_option_descrips = """
        -w  treat duplicates in barcode_file as warning, not error
        -fh first line of barcode file considered a header (default: auto detect header)
        -nh first line of barcode file is not a header (default: auto detect header)
        -info cols show column settings in barcode file and values for the first line
        -info all|fwd|rev|both|primer display barcode index or primer info, including edit distances

        -n <num_seqs> number of sequences to read from file (ex: -n 100)
        -n <first_seq>,<num_seqs> (ex: -n 102,3)\n"""

    extra = secondary_option_descrips if show_all_descrips else ""
    extra_fmts = secondary_output_formats if show_all_descrips else ""

    usage = """
    Usage: minibar.py barcode_file sequence_file [-pct <pct> | -e <int> -E <int>] [-l <int>]
                                                 [-F [-P <prefix>]] [-M 1|2|3]
                                                 [-S | -T | -C | -CC | -D]
                                                 [-cols <int_list>] [-info cols|fwd|rev|primer]
                                                 [-w] [-fh | -nh] [-n <num_seqs> | -n <first_seq>,<num_seqs>]

        Identify MinION sequence by dual barcode indexes and primers.
        The sequence file can be in Fasta or Fastq format, gzipped or plain text.
        Sample ID is placed at end of header comment with match hit info before it.
        ({})

        Example: ./minibar.py -C -F Demultiplex.txt example.fq

        -h display this with all option's descriptions     -v displays version
        -p <pct> percentage match (.75)
        -e <int> barcode edit distance value, overrides -p (4)
        -E <int> primer edit distance value (11)
        -l <int> length to search for index and primer at start and end of sequence (80)

        -F create individual sample files for sequences with -S or -C output (default: False)
        -f write to stdout instead of creating files (default: True)
        -P <str> if -F, <str> is prefix for individual files, followed by sample ID. (default: sample_)

        -M 1|2|3 Method to identify sample types using the barcodes (default: 3)
                 1 requires approximate match of barcode and primer at sequence start, this
                   and barcodes matched at the other end are used to identify sample IDs
                 2 finds matched barcodes on both ends of sequence, identifies pairs that match a sample ID
                 3 uses Method 1 and if it does not succeed, uses Method 2
                (0 identify sample ID even if only one barcode found -- ambiguous if barcodes not unique)

        -S outputs sequence record in fasta or fastq format of input (default output)
        -T trims barcode and primer from each end of the sequence, then outputs record
        -C similar to S but uses upper/lower case to show found barcode indexes and primers
        -CC also colors found barcode blue, primer green if found, primer red otherwise
        {}
        -cols <int_list> column position in barcode_file for: sample, fwd index, fwd primer, rev index, rev primer
                 (default: 1,2,3,4,5 if 5 cols; 1,3,4,5,6 if 6 cols; 1,3,4,6,7 if 7 cols; 1,3,5,8,10 if 10 or more cols)
        {}""".format(version(), extra_fmts, extra)

    print(usage, file=sys.stderr)
    sys.exit(3)


def getoptions(argv):
    def is_int(i):
        try:
            int(i)
            return True
        except:
            return False

    def is_float(f):
        try:
            float(f)
            return True
        except:
            return False

    def num_err(arg, val):
        error('Expected numbers, got \'{} {}\''.format(arg, val))

    def process_cols_list(cols):
        num_cols = 5 if opts.dual_indexes else 3
        lst = cols.split(',')
        if len(lst) != num_cols:
            error('-cols value requires {} comma delimited positions.'.format(num_cols))
        for i in range(num_cols):
            if not is_int(lst[i]) or int(lst[i])<1:
                error('-cols value requires {} comma delimited numbers greater than 1.'.format(num_cols))

        ilst = [int(i) for i in lst]

        opts.sampleix = ilst[0]-1
        opts.fwix = ilst[1]-1
        opts.fwprm = ilst[2]-1
        if num_cols >= 5:
            opts.rvix = ilst[3]-1
            opts.rvprm = ilst[4]-1
        else:
            opts.rvix = -1
            opts.rvprm = -1
        opts.col_set_on_cmdline = True

    min_options_required = 2
    num_args = len(argv)
    if num_args < min_options_required:
        usage()

    class opts:
        primerfile = ''
        sequencefile = ''

        dual_indexes = True
        search_len = 80
        num_seqs = -1  # -1 means all sequences
        start_seq = 0  # 0 means start with first seq
        index_edit_distance = .75  # if < 0 it is a percentage to calculate
        primer_edit_distance = -1  # -1 means calculate based on primer length
        search_method = 3  # 2: means weak match of pairs, 1: strong ix/prm match at seq begin, 3: try 1, then 2

        output_to_files = False  # True means write sequences to individual sample ID files
        output_file_prefix = 'sample_'
        isfastq = True  # reset when file sequence actually read
        output = 1  # 4 trim, 3 is diagnostic, 2 means use case to show index & primers in sequence, else as is
        output_letter = 'S'  # letter that corresponds to S
        show_color = False  # in output 2 show color along with change in case (-F turns off)

        # FWIX = 2; RVIX = 7; FWPRM = 4; RVPRM = 9; SAMPLE = 0
        sampleix = 0
        fwix     = 2
        fwprm    = 4
        rvix     = 7
        rvprm    = 9
        error_on_duplicate_samples = True
        col_set_on_cmdline = False  # if using defaults do settings based on num cols in file
        first_line_header = -1  # -1 auto detect, 0 don't treat as header, otherwise it's a header
        barcode_file_info = ''  # valid values start with either 'f' 'r' or 'p'

    unrecog = ''
    ix_arg = 0
    while ix_arg < (num_args - 1):
        ix_arg += 1
        arg = argv[ix_arg]
        if arg == '-h' or arg == '--help':
            usage(True)
        elif arg == '-v' or arg == '--version':
            print(version())
            sys.exit(0)

        if arg[0] != '-':  # first 2 args without '-' are primer file and sequence file name
            if opts.primerfile == "":
                opts.primerfile = arg
            elif opts.sequencefile == "":
                opts.sequencefile = arg
            else:
                unrecog = arg
        elif len(arg) > 1:
            op = arg[1]
            if op in 'SCDT':  # single letter options
                opts.output = {'S':1, 'C':2, 'D':3, 'T':4}.get(op,1)
                opts.output_letter = op
                if op == 'C' and len(arg)>2: opts.show_color = True
            elif arg == '-nh' or arg == '-fh':  # args to force 1st line as header or force as non-header
                opts.first_line_header = 0 if arg == '-nh' else '1'
            elif op == 'F':  # trigger writing each sample to its own file
                opts.output_to_files = True
            elif op == 'f':  # trigger writing each sample to its own file
                opts.output_to_files = False
            elif op == 'w':  # treat duplicate samples as warning not error
                opts.error_on_duplicate_samples = False
            elif op in 'pkeKElMncPi':  # these have a value after option
                ix_arg += 1
                if ix_arg > (num_args - 1):
                    error('expected value after {}'.format(arg))
                val = argv[ix_arg]
                have_int = is_int(val); have_float = is_float(val)
                if op == 'n':  # -n has 2 formats -n int and -n int,int
                    if ',' in val:
                        n1, n2 = val.split(',')
                        if not is_int(n1) or not is_int(n2): num_err(arg, val)
                        opts.start_seq = int(n1)
                        opts.num_seqs  = int(n2)
                    else:
                        if not have_int: num_err(arg, val)
                        opts.num_seqs = int(val)
                elif op == 'c':
                    if arg == '-cols':
                        process_cols_list(val)
                    else:
                        unrecog = arg
                elif op == 'i':
                    if arg == '-info':
                        opts.barcode_file_info = val
                    else:
                        unrecog = arg

                if (op in 'lkeKEM' and not have_int) or (op == 'p' and not have_float):
                    num_err(arg, val)

                if op == 'l':
                    opts.search_len = int(val)
                elif op == 'p':
                    opts.index_edit_distance = float(val)
                elif op == 'k' or op == 'e':
                    opts.index_edit_distance = int(val)
                elif op == 'K' or op == 'E':
                    opts.primer_edit_distance = int(val)
                elif op == 'P':
                    opts.output_file_prefix = val
                elif op == 'M':
                    opts.search_method = int(val) if int(val) in [1,2,3,0] else 1
            else:
                unrecog = arg
                break

    if opts.output_to_files:
        opts.show_color = False
    if unrecog != '':
        error('"{}" is an unrecognized option'.format(unrecog))
    elif opts.primerfile == "":
        error('missing barcode file')
    elif opts.barcode_file_info != '':
        display_barcode_file_inf(opts)
        sys.exit(0)
    elif opts.sequencefile == "":
        error('missing sequence file')

    return opts


def main(argv):  # pass in argv in case we want to use as a module and pass in arg list
    try:
        global rev_comp # use translate table appropriate for the major version of Python
        python_major = sys.version_info[0]
        rev_comp = rev_comp_py2 if python_major < 3 else rev_comp_py3

        options = getoptions(argv)
        minibar(options)

    except (KeyboardInterrupt, SystemExit):
        sys.exit(0)
    except IOError as e:
        if e.errno == errno.EPIPE:
            # stdout is closed, no point in continuing
            # attempt to close stdout, stderr explicitly to prevent cleanup problems:
            try:
                sys.stdout.close()
            except IOError:
                pass
            try:
                sys.stderr.close()
            except IOError:
                pass
        elif e.errno == errno.ENOENT:
            error('No file named \'{}\''.format(e.filename))
        else:  # show traceback when IO error other than pipe closed or file not found
            raise


if __name__ == "__main__":
    main(sys.argv)
