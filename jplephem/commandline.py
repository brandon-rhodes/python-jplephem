"""The `python -m jplephem` command line."""

from __future__ import print_function

import argparse
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
    p.set_defaults(func=print_daf_segments)
    p.add_argument('path', help='Path to a SPICE file')

    p = subparsers.add_parser(
        'spk',
        help="List the segments in an SPK file",
    )
    p.set_defaults(func=print_spk_segments)
    p.add_argument('path', help='Path to a .bsp SPICE kernel file')

    args = parser.parse_args(args)
    func = getattr(args, 'func', None)
    if func is None:
        parser.print_help()
        return 2

    func(args)
    return 0

def print_daf_segments(args):
    with open(args.path, 'rb') as f:
        d = DAF(f)
        for i, (name, values) in enumerate(d.summaries()):
            print('%2d' % (i+1),
                  name.decode('latin-1'),
                  ' '.join(repr(v) for v in values))

def print_spk_segments(args):
    with open(args.path, 'rb') as f:
        print(str(SPK(DAF(f))))
