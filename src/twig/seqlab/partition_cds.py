#!/usr/bin/env python

"""...

awk 'NR>1 && /^>/ {exit} !/^>/ {L+=length($0)} END{
  print "GTR+G, gene_pos1 = 1-"L"\\3"
  print "GTR+G, gene_pos2 = 2-"L"\\3"
  print "GTR+G, gene_pos3 = 3-"L"\\3"
}' "${PRE}.final.nt.fa" > "${PRE}.final.nt.fa.partitions.txt"
"""


def run_partition_cds(args):
    """Runs diamond blast with args parsed from sys or CLI."""
    if not args.input.exists():
        raise IOError(f"{args.input} does not exist")
    with args.input.open('rt') as hin:
        lines = hin.readlines()
        length = max(len(i.strip()) for i in lines if not i.startswith(">"))
    if length % 3:
        raise ValueError("CDS length is not a multiple of 3")

    # require nmodels to match n...
    if len(args.model) == 1:
        args.model = args.model * args.number
    if args.number != len(args.model):
        raise ValueError("--number should match length of --model")

    if args.number == 2:
        print(f"{args.model[0]}, gene_pos1_2 = 1-{length}\\3, 2-{length}\\3")
        print(f"{args.model[1]}, gene_pos3 = 1-{length}\\3")

    else:
        print(f"{args.model[0]}, gene_pos1 = 1-{length}\\3")
        print(f"{args.model[1]}, gene_pos2 = 2-{length}\\3")
        print(f"{args.model[2]}, gene_pos3 = 3-{length}\\3")


def main():
    from twig.cli.subcommands import get_parser_partition_cds
    parser = get_parser_partition_cds()
    args = parser.parse_args()
    run_partition_cds(args)


if __name__ == "__main__":
    main()
