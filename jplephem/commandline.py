"""The `python -m jplephem` command line."""

from __future__ import print_function

import argparse
import sys
from .daf import DAF
from .excerpter import RemoteFile, write_excerpt
from .spk import SPK

def main(args):
    parser = argparse.ArgumentParser(
        prog='python -m jplephem',
        description='Describe an SPK kernel',
    )
    subparsers = parser.add_subparsers()

    p = subparsers.add_parser(
        'comment',
        help="Print a file's comment blocks",
    )
    p.set_defaults(func=comment)
    p.add_argument('path', help='Path to a SPICE file')

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
    p.add_argument('--targets', help='Comma-separated targets to include')
    p.add_argument('start_date', help='Start date yyyy/mm/dd', type=parse_date)
    p.add_argument('end_date', help='End date yyyy/mm/dd', type=parse_date)
    p.add_argument('path_or_url', help='Local filename or remote URL')
    p.add_argument('output_path', help='Output file to create')

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
    if lines and not lines[-1].endswith('\n'):
        lines.append('')
    return '\n'.join(lines)

def comment(args):
    with open(args.path, 'rb') as f:
        d = DAF(f)
        yield d.comments()

def daf_segments(args):
    with open(args.path, 'rb') as f:
        d = DAF(f)
        for i, (name, values) in enumerate(d.summaries()):
            yield '{:2d} {} {}'.format(i + 1, name.decode('latin-1'),
                                       ' '.join(repr(v) for v in values))

def excerpt(args):
    if args.path_or_url.startswith(('http://', 'https://')):
        url = args.path_or_url
        f = RemoteFile(url)
    else:
        path = args.path_or_url
        f = open(path, 'rb')

    with f:
        spk = SPK(DAF(f))
        summaries = spk.daf.summaries()

        if args.targets:
            desired_targets = set(args.targets.split(','))
            summaries = [
                summary for summary, segment in zip(summaries, spk.segments)
                if str(segment.target) in desired_targets
            ]

        with open(args.output_path, 'w+b') as output_file:
            write_excerpt(spk, output_file, args.start_date, args.end_date,
                          summaries)

    return ()

def spk_segments(args):
    with open(args.path, 'rb') as f:
        yield str(SPK(DAF(f)))

def parse_date(s):
    try:
        fields = [int(f) for f in s.split('/')]
    except ValueError:
        fields = []
    if len(fields) < 1 or len(fields) > 3:
        E = argparse.ArgumentTypeError
        raise E('specify each date as YYYY or YYYY/MM or YYYY/MM/DD')
    return julian_day(*fields)

def filter_segments(spk, segment_names):
    segment_names = set(segment_names)
    spk.segments = [
        segment for segment in spk.segments
        if str(segment.target) in segment_names
    ]

def julian_day(year, month=1, day=1):
    """Given a proleptic Gregorian calendar date, return a Julian day int."""
    janfeb = month < 3
    return (day
            + 1461 * (year + 4800 - janfeb) // 4
            + 367 * (month - 2 + janfeb * 12) // 12
            - 3 * ((year + 4900 - janfeb) // 100) // 4
            - 32075)
