"""The `python -m jplephem` command line."""

from __future__ import print_function

import argparse
import sys
from .daf import DAF
from .spk import SPK

def main(args):
    parser = argparse.ArgumentParser(
        prog='python -m jplephem',
        description='Describe an SPK kernel',
    )
    subparsers = parser.add_subparsers()

    p = subparsers.add_parser(
        'daf',
        help="List a file's raw segment descriptors",
    )
    p.set_defaults(func=daf_segments)
    p.add_argument('path', help='Path to a SPICE file')

    p = subparsers.add_parser(
        'excerpt',
        help="Create an SPK covering a narrower range of dates",
    )
    p.set_defaults(func=excerpt)
    p.add_argument('start_date', help='Start date yyyy/mm/dd')
    p.add_argument('end_date', help='End date yyyy/mm/dd')
    p.add_argument('file', help='Local filename or remote URL')

    p = subparsers.add_parser(
        'spk',
        help="List the segments in an SPK file",
    )
    p.set_defaults(func=spk_segments)
    p.add_argument('path', help='Path to a .bsp SPICE kernel file')

    args = parser.parse_args(args)
    func = getattr(args, 'func', None)
    if func is None:
        parser.print_help()
        sys.exit(2)

    lines = list(func(args))
    lines.append('')
    return '\n'.join(lines)

def daf_segments(args):
    with open(args.path, 'rb') as f:
        d = DAF(f)
        for i, (name, values) in enumerate(d.summaries()):
            yield '{:2d} {} {}'.format(i + 1, name.decode('latin-1'),
                                       ' '.join(repr(v) for v in values))

def excerpt(args):
    pass

def spk_segments(args):
    with open(args.path, 'rb') as f:
        yield str(SPK(DAF(f)))
