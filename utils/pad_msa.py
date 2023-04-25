
import argparse

# def pad_a3m(in_fname, out_fname, n_g):
#     # add poly g to multimer msa
#     infile = open(in_fname, 'r')
#     outfile = open(out_fname, 'w')

#     count = 0
#     space = "-" * n_g
#     poly_g = "G" * n_g
#     eof, add_space, add_poly_g = False, False, False

#     while not eof:
#         count += 1
#         line = infile.readline()
#         if len(line) == 0:
#             out_line = line
#             eof = True

#         elif count == 1:
#             assert(line[0] == '#')
#             segments = line[1:].split()
#             len1, len2 = segments[0].split(',')
#             len1 = int(len1)
#             len2 = int(len2)
#             out_line = f"#{len1+n_g},{len2}   {segments[1]}\t"

#         elif count == 2:
#             assert(line[0] == '>')
#             id1, id2 = line[1:].split()
#             add_poly_g = True
#             out_line = line

#         elif line[0] == '>':
#             add_space = True
#             out_line = line

#         else:
#             assert(add_space ^ add_poly_g)
#             if add_poly_g:
#                 out_line = line[:len1] + poly_g + line[len1+1:]
#                 add_poly_g = False
#             elif add_space:
#                 out_line = line[:len1] + space + line[len1+1:]
#                 add_space = False

#         outfile.writelines(out_line)

#     infile.close()
#     outfile.close()

def pad_a3m(in_fname, out_fname, n_g):
    infile = open(in_fname, 'r')
    outfile = open(out_fname, 'w')

    count = 0
    space = "-" * n_g
    poly_g = "G" * n_g
    eof, add_space, replace_full_seq = False, False, False

    # read 1st line
    line = infile.readline()
    assert(line[0] == '#')
    segments = line[1:].split()
    len1, len2 = segments[0].split(',')
    len1 = int(len1)
    len2 = int(len2)
    suffix = segments[1].split(',')[0]
    out_line = f"#{len1+n_g+len2}\t{suffix}\n" # keep only one chain
    outfile.writelines(out_line)

    # read 2nd line (dont write to output)
    line = infile.readline()
    assert(line[0] == '>')
    id1, id2 = line[1:].split()

    # read 3rd line (dont write to output)
    line = infile.readline()
    full_seq = line[:len1] + poly_g + line[len1+1:]

    # read rest lines
    while not eof:
        count += 1
        line = infile.readline()
        if len(line) == 0:
            out_line = line
            eof = True

        elif line[0] == '>':
            if line[1:4] == id1 or line[1:4] == id2:
                replace_full_seq = True
                out_line = '>' + id1 + '\n'
            else:
                add_space = True
                out_line = line
        else:
            assert(add_space ^ replace_full_seq)
            if replace_full_seq:
                out_line = full_seq
                replace_full_seq = False
            elif add_space:
                out_line = line[:len1] + space + line[len1:]
                add_space = False

        if out_line is not None:
            outfile.writelines(out_line)

    infile.close()
    outfile.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-g", type=int)
    parser.add_argument("--in-fname", type=str)
    parser.add_argument("--out-fname", type=str)

    args = parser.parse_args()
    pad_a3m(args.in_fname, args.out_fname, args.n_g)
