"""Extract data for a specific date range from an SPK file."""

from sys import stderr
try:
    from urllib.request import URLopener
except:
    from urllib import URLopener

from numpy import copy
from .daf import DAF
from .spk import S_PER_DAY, T0

clip_lower = max
clip_upper = min

def _seconds(jd):
    """Convert a Julian Date to a number of seconds since J2000."""
    return (jd - T0) * S_PER_DAY

def write_excerpt(input_spk, output_file, start_jd, end_jd, summaries):
    start_seconds = _seconds(start_jd)
    end_seconds = _seconds(end_jd)
    old = input_spk.daf

    # Copy the file record and the comments verbatim.
    f = output_file
    f.seek(0)
    f.truncate()
    for n in range(1, old.fward):
        data = old.read_record(n)
        f.write(data)

    # Start an initial summary and name block.
    summary_data = b'\0' * 1024
    name_data = b' ' * 1024
    f.write(summary_data)
    f.write(name_data)

    d = DAF(f)
    d.fward = d.bward = old.fward
    d.free = (d.fward + 1) * (1024 // 8) + 1
    d.write_file_record()

    # Copy over an excerpt of each array.
    for name, values in summaries:
        start, end = values[-2], values[-1]
        init, intlen, rsize, n = old.read_array(end - 3, end)
        rsize = int(rsize)

        i = int(clip_lower(0, (start_seconds - init) // intlen))
        j = int(clip_upper(n, (end_seconds - init) // intlen + 1))
        init = init + i * intlen
        n = j - i

        extra = 4     # enough room to rebuild [init intlen rsize n]
        excerpt = copy(old.read_array(
            start + rsize * i,
            start + rsize * j + extra - 1,
        ))
        excerpt[-4:] = (init, intlen, rsize, n)
        values = (init, init + n * intlen) + values[2:]
        d.add_array(b'X' + name[1:], values, excerpt)

class RemoteFile(object):
    def __init__(self, url):
        self.opener = URLopener()
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
        h = 'Range', 'bytes={}-{}'.format(start, end)
        stderr.write('Fetching {} {}\n'.format(self.filename, h[1]))
        self.opener.addheaders.append(h)
        data = self.opener.open(self.url).read()
        return data

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
