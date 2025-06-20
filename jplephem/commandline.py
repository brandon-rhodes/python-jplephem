"""The `python -m jplephem` command line."""

from __future__ import print_function

import argparse
import sys
from .calendar import compute_calendar_date, compute_julian_date
from .daf import DAF
from .excerpter import RemoteFile, write_excerpt
from .spk import S_PER_DAY, SPK, T0

_DAY = 86400.0

def _jd(seconds):
    """Convert a number of seconds since J2000 to a Julian Date."""
    return T0 + seconds / S_PER_DAY

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
    p.add_argument('-v', '--verbose', action='store_true')

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
    for string, jd in args.start_date, args.end_date:
        yield 'Date {:10} = JD {}'.format(string, jd)

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
            write_excerpt(spk, output_file, args.start_date[1],
                          args.end_date[1], summaries)

    yield '\n{!r} written successfully with the following contents\n'.format(
        args.output_path)

    with open(args.output_path, 'rb') as f:
        yield str(SPK(DAF(f)))

def spk_segments(args):
    with open(args.path, 'rb') as f:
        spk = SPK(DAF(f))

        # Snag the first line from the normal str().

        output = str(spk)
        yield output.split('\n', 1)[0]

        # But produce the rest of the lines ourselves, so we can
        # optionally honor '-v' by providing more information.

        for s in spk.segments:
            yield str(s)
            if not args.verbose:
                continue
            for line in _describe_segment_details(s):
                yield line

def _describe_segment_details(s):
    if s.data_type not in (2, 3):
        return

    init, intlen, coefficients = s._data
    degree, dimensions, record_count = coefficients.shape
    days_per = intlen / _DAY

    plural = '' if record_count == 1 else 's'
    each = '' if record_count == 1 else ' each'
    yield '   {} polynomial{} covering {} days{}'.format(
        record_count, plural, days_per, each,
    )
    yield '      x {} coefficients per polynomial'.format(degree)
    yield '      x {} coordinates'.format(dimensions)
    yield '      = {} double precision floats'.format(coefficients.size)

    polynomial_start = init
    polynomial_end = init + intlen * record_count
    if s.start_second == polynomial_start:
        yield '   Polynomial start date matches segment start date'
    else:
        days = (s.start_second - polynomial_start) / _DAY
        jd = _jd(polynomial_start)
        y, m, d = compute_calendar_date(int(jd + 0.5))
        yield (
            '   First polynomial starts {:.1f} days earlier'
             ' than segment start date, on {}-{:02}-{:02}'
             .format(days, y, m, d)
        )
    if s.end_second == polynomial_end:
        yield '   Polynomial end date matches segment end date'
    else:
        days = (polynomial_end - s.end_second) / _DAY
        jd = _jd(polynomial_end)
        y, m, d = compute_calendar_date(int(jd + 0.5))
        yield (
            '   Final polynomial ends {:.1f} days later'
            ' than segment end date, on {}-{:02}-{:02}'
            .format(days, y, m, d)
        )

    yield ''

def parse_date(s):
    try:
        fields = [int(f) for f in s.split('/')]
    except ValueError:
        fields = []
    if len(fields) < 1 or len(fields) > 3:
        E = argparse.ArgumentTypeError
        raise E('specify each date as YYYY or YYYY/MM or YYYY/MM/DD')
    jd = compute_julian_date(*fields)
    return s, jd
