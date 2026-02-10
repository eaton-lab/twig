#!/usr/bin/env python

"""Run MACSE codon-aware CDS alignment."""

from loguru import logger

from .align_macse import run_align_macse as run_align_cds


def main() -> None:
    from twig.cli.subcommands import get_parser_align_cds

    parser = get_parser_align_cds()
    args = parser.parse_args()
    run_align_cds(args)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("interrupted by user")
    except Exception as exc:
        logger.error(exc)
