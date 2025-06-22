"""Extract data for a specific date range from an SPK file."""

from sys import stderr
try:
    from urllib.request import Request, urlopen
except ImportError:
    from urllib2 import Request, urlopen

from numpy import copy
from . import __version__ as jplephem_version
from .calendar import compute_calendar_date
from .daf import DAF, K
from .spk import S_PER_DAY, T0

clip_lower = max
clip_upper = min

_PREFACE = """\
;
; This is an ephemeris excerpt created by jplephem {}, which was
; asked to narrow the ephemeris to Julian dates {:.1f} - {:.1f}
; (proleptic Gregorian dates {}-{:02}-{:02} through {}-{:02}-{:02}).
;
; Here is the comments area from the original ephemeris file:
; ----------------------------------------------------------------------
"""

def _seconds(jd):
    """Convert a Julian Date to a number of seconds since J2000."""
    return (jd - T0) * S_PER_DAY

def write_excerpt(input_spk, output_file, start_jd, end_jd, summaries):
    start_seconds = _seconds(start_jd)
    end_seconds = _seconds(end_jd)
    old = input_spk.daf

    # Supplement the comment text.
    y1, m1, d1 = compute_calendar_date(int(start_jd + 0.5))
    y2, m2, d2 = compute_calendar_date(int(end_jd + 0.5))
    preface = _PREFACE.format(
        jplephem_version, start_jd, end_jd, y1,m1,d1, y2,m2,d2,
    )
    comment = preface + old.comments()

    # Build new comment blocks (which have 1000 text characters each).
    data = comment.encode('ascii').replace(b'\n', b'\0') + b'\004'
    blocks = [data[i : i + 1000] for i in range(0, len(data), 1000)]
    comment_data = b''.join([
        block + b' ' * (K - len(block))
        for block in blocks
    ])

    # Start the new DAF file with:
    # 1. The verbatim first record from the original file.
    # 2. The new comment.
    # 3. An empty summary block.
    # 4. An empty name block.

    f = output_file
    f.seek(0)
    f.truncate()

    summary_data = b'\0' * 1024
    name_data = b' ' * 1024

    f.write(old.read_record(1))
    f.write(comment_data)
    f.write(summary_data)
    f.write(name_data)

    # There are now enough blocks to start treating the file as a DAF!
    # Set the initial block number indexes.
    f.seek(0)
    d = DAF(f)
    d.fward = d.bward = 2 + len(comment_data) // K
    d.free = (d.fward + 1) * (1024 // 8) + 1
    d.write_file_record()

    # Copy over an excerpt of each array.
    for name, values in summaries:
        start, end = values[-2], values[-1]
        init, intlen, rsize, n = old.read_array(end - 3, end)
        rsize = int(rsize)

        i = int(clip(0, n, (start_seconds - init) // intlen))
        j = int(clip(0, n, (end_seconds - init) // intlen + 1))
        if i == j:
            continue  # Segment has no overlap with user's dates.

        init = init + i * intlen
        n = j - i

        extra = 4     # enough room to rebuild [init intlen rsize n]
        excerpt = copy(old.read_array(
            start + rsize * i,
            start + rsize * j + extra - 1,
        ))
        excerpt[-4:] = (init, intlen, rsize, n)

        # Even though the polynomials we selected probably cover a wider
        # range of dates, let's only claim that each segment covers the
        # range `start_seconds .. end_seconds` that the user requested,
        # to avoid confusing them.  The `de442s.bsp` ephemeris shows
        # that this is also the practice at NASA itself.

        values = (start_seconds, end_seconds) + values[2:]
        d.add_array(name, values, excerpt)

def clip(lower, upper, n):
    return clip_lower(lower, clip_upper(upper, n))

class RemoteFile(object):
    def __init__(self, url):
        self.url = url
        self.filename = url.rstrip('/').rsplit('/', 1)[-1]
        self.offset = 0

    def seek(self, offset, whence=0):
        assert whence == 0
        self.offset = offset

    def read(self, size):
        start = self.offset
        end = start + size - 1
        assert end > start
        byte_range = 'bytes={}-{}'.format(start, end)
        stderr.write('Fetching {} bytes from {} using Range: {}\n'
                     .format(size, self.filename, byte_range))
        request = Request(self.url, headers={'Range': byte_range})
        data = urlopen(request).read()
        assert len(data) == size, (
            'asked for "Range: {}" which is {} bytes, but got {} bytes back'
            .format(byte_range, size, len(data))
        )
        self.offset += size
        return data

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
