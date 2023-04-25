
import argparse

def pad_a3m(in_fname, out_fname, n_g):
    infile = open(in_fname, 'r')
    outfile = open(out_fname, 'w')

    count = 0
    space = "-" * n_g
    poly_g = "G" * n_g
    eof, add_space, add_poly_g = False, False, False

    while not eof:
        count += 1
        line = infile.readline()
        if len(line) == 0:
            out_line = line
            eof = True

        elif count == 1:
            assert(line[0] == '#')
            segments = line[1:].split()
            len1, len2 = segments[0].split(',')
            len1 = int(len1)
            len2 = int(len2)
            out_line = f"#{len1+n_g},{len2}   {segments[1]}\n"

        elif count == 2:
            assert(line[0] == '>')
            id1, id2 = line[1:].split()
            add_poly_g = True
            out_line = line

        elif line[0] == '>':
            add_space = True
            out_line = line

        else:
            assert(add_space ^ add_poly_g)
            if add_poly_g:
                out_line = line[:len1] + poly_g + line[len1+1:]
                add_poly_g = False
            elif add_space:
                out_line = line[:len1] + space + line[len1+1:]
                add_space = False

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
