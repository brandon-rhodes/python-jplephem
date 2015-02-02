"""Interpret a binary DAF file.

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html

"""
import mmap
import struct
from pprint import pprint

BFF = 'BIG-IEEE', 'LTL-IEEE', 'VAX-GFLT', 'VAX-DFLT'      # Binary file format
FTPSTR = 'FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP'  # FTP test string
RECORD_LENGTH = 1024

class DAF(object):
    def __init__(self, open_file):
        self.map = mmap.mmap(open_file.fileno(), 0, access=mmap.ACCESS_READ)
        self.read_file_record()

    def read_file_record(self):
        map = self.map
        self.locfmt = map[88:96]

        if self.locfmt == b'BIG-IEEE':
            self.endian = '>'
        elif self.locfmt == b'LTL-IEEE':
            self.endian = '<'
        else:
            raise ValueError('unrecognized format: {0!r}'.format(self.locfmt))

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

        if map[500:1000].strip(b'\0') != FTPSTR:
            raise ValueError('the file has been damaged')

        n = self.fward - 1
        i = RECORD_LENGTH * n

        summary_record = map[i:i + RECORD_LENGTH]
        name_record = map[i + RECORD_LENGTH:i + 2 * RECORD_LENGTH]

        next_summary, previous_summary, n_summaries = struct.unpack(
            self.endian + 'ddd', summary_record[:24])

        n_summaries = int(n_summaries)
        print n_summaries

        step = self.summary_step
        length = self.summary_length

        for i in range(0, n_summaries * step, step):
            print name_record[i:i+step],
            j = 24 + i
            print struct.unpack(self.summary_format, summary_record[j:j+length])

        pprint(vars(self))


def main():
    with open('jup310.tmp', 'rb') as f:
        daf = DAF(f)

    # print repr(b[128:128+32])
    # print b.index('LTL-IEEE')

main()
