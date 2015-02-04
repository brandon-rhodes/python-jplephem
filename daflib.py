"""Interpret a binary DAF file.

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html
http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

"""
import mmap
import numpy as np
import struct
import sys
from collections import namedtuple
from jplephem.ephem import Ephemeris

BFF = b'BIG-IEEE', b'LTL-IEEE', b'VAX-GFLT', b'VAX-DFLT'   # Binary file format
FTPSTR = b'FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP'  # FTP test string
RECORD_LENGTH = 1024
S_PER_DAY = 86400.0
T0 = 2451545.0

Summary = namedtuple('Summary', 'source start_second stop_second target'
                     ' center frame data_type start_index stop_index')

class DAF(object):
    """Access to NASA SPICE Double Precision Array Files (DAF)."""

    def __init__(self, open_file):
        self.map = mmap.mmap(open_file.fileno(), 0, access=mmap.ACCESS_READ)
        if sys.version_info > (3,):
            self.map = memoryview(self.map)
        self.read_file_record()

    def read_file_record(self):
        map = self.map
        self.locfmt = map[88:96]

        if self.locfmt == b'BIG-IEEE':
            self.endian = '>'
        elif self.locfmt == b'LTL-IEEE':
            self.endian = '<'
        else:
            raise ValueError('unrecognized format {0!r}'.format(self.locfmt))

        (locidw, self.nd, self.ni, locifn, self.fward, self.bward,
         self.free) = struct.unpack(self.endian + '8sII60sIII', map[:88])

        self.ss = self.nd + (self.ni + 1) // 2

        self.summary_step = 8 * self.ss
        self.summary_format = self.endian + 'd' * self.nd + 'i' * self.ni
        self.summary_length = struct.calcsize(self.summary_format)

        self.locidw = locidw.upper().rstrip()
        self.locifn = locifn.upper().rstrip()

        if self.locidw != b'DAF/SPK':
            raise ValueError(
                'the first bytes do not identify this as a DAF/SPK file')

        if bytes(map[500:1000]).strip(b'\0') != FTPSTR:
            raise ValueError('the file has been damaged')

    def record(self, n):
        """Return record `n` as bytes, where the first record is record 1."""
        start = RECORD_LENGTH * (n - 1)
        return self.map[start:start + RECORD_LENGTH]

    def summaries(self):
        """Yield a (name, values) tuple for each summary in the file."""

        record_number = self.fward
        length = self.summary_length
        step = self.summary_step

        while record_number:
            summary_record = self.record(record_number)
            name_record = self.record(record_number + 1)

            next_number, previous_number, n_summaries = struct.unpack(
                self.endian + 'ddd', summary_record[:24])

            for i in range(0, int(n_summaries) * step, step):
                source = bytes(name_record[i:i+step]).strip()
                j = 24 + i
                data = summary_record[j:j+length]
                values = struct.unpack(self.summary_format, data)
                yield Summary(source, *values)

            record_number = int(next_number)

    def bytes(self, start, stop):
        """Return data from double-precision word start to stop, inclusive."""
        return self.map[8 * start - 8 : 8 * stop]

    def __getitem__(self, index):
        if isinstance(index, slice):
            start = index.start
            stop = index.stop
            format = self.endian + 'd' * (1 + stop - start)
            print(8 * start)
            return struct.unpack(format, self.map[8 * start - 8:8 * stop])
        return struct.unpack(self.endian + 'd', self.map[
            8 * index - 8, 8 * index])


class SPK(Ephemeris):
    """A JPL SPK-format ephemeris that computes positions and velocities."""

    def __init__(self, path):
        self.daf = DAF(open(path))
        self.sets = {}
        self.summaries = {s.target: s for s in self.daf.summaries()}

    def load(self, target):
        s = self.sets.get(target)
        if s is None:
            summary = self.summaries[target]
            if summary.data_type == 2:
                component_count = 3
            elif summary.data_type == 3:
                component_count = 6
            else:
                raise ValueError('only SPK data types 2 and 3 are supported')
            stop = summary.stop_index
            init, intlen, rsize, n = self.daf[stop-3:stop]
            coefficient_count = (rsize - 2) // component_count
            print(summary)
            # TODO: use intlen directly to create days_per_set
            self.jalpha = T0 + init / S_PER_DAY
            self.jomega = self.jalpha + intlen * n / S_PER_DAY
            print('omega:', init + intlen * n)
            print(init, intlen, rsize, n)
            data = self.daf.bytes(summary.start_index, stop)
            s = np.ndarray((n, rsize), self.daf.endian + 'd', data)
            s = s[:,2:]
            s.shape = (n, component_count, coefficient_count)
        return s


def main2():
    T0 = 2451545.0
    s = SPK('jup310.bsp')
    p = s.position(3, T0)
    print(p)
    p = s.position(502, T0)
    print(p)

def main():
    with open('jup310.bsp', 'rb') as f:
        daf = DAF(f)

    for name, values in daf.summaries():
        print(name, values)
        (initial_epoch, final_epoch, target_code, center_code, frame_code,
         data_type, start, end) = values
        break

    n = 7208500
    init, intlen, rsize, n = daf[n-3:n+1]
    rsize = int(rsize)
    n = int(n)
    # print(rsize)
    # print(rsize * n)
    # print(7208500 + 1 - 897 - 4)
    coefficient_count = (rsize - 2) // 6  # -2 for mid, radius
    # print(coefficient_count)

    mid = numpy.ndarray(
        (8,),
        daf.endian + 'd',
        daf.bytes(start, end + 1),
        strides=rsize * 8,
    )
    print(mid)

    radius = numpy.ndarray(
        (8,),
        daf.endian + 'd',
        daf.bytes(start, end + 1),
        offset=8,
        strides=rsize * 8,
    )
    print(radius)

    radius = numpy.ndarray(
        (8,),
        daf.endian + 'd',
        daf.bytes(start, end + 1),
        offset=8,
        strides=rsize * 8,
    )
    print(radius)

    a = numpy.ndarray(
        (n, rsize),
        daf.endian + 'd',
        daf.bytes(start, end + 1),
    )
    print('-------')
    print(a[0])
    print('-------')
    # print(a[1])
    b = a[:,0]
    c = a[:,1]
    d = a[:,2::coefficient_count]
    e = a[:,3::coefficient_count]
    f = a[:,2+coefficient_count-1::coefficient_count]
    print(b)
    print(c)
    print(d[0])
    print(e[0])
    print(f[0])
    #b = a.resize
    g = a[:,2:]
    g.shape = (n, 6, coefficient_count)
    #h = g.reshape(
    print('==========')
    for i in range(0, coefficient_count):
        print(g[0,:,i])
    # print repr(b[128:128+32])
    # print b.index('LTL-IEEE')

#main()
main2()
