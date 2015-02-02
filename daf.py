"""Interpret a binary DAF file.

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html

"""
import struct

FTPSTR = 'FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP'  # FTP test string
BFF = 'BIG-IEEE', 'LTL-IEEE', 'VAX-GFLT', 'VAX-DFLT'      # Binary file format

class DAF(object):
    def __init__(self, open_file):
        self.f = open_file
        self.read_file_record()

    def read_file_record(self):
        f = self.f
        f.seek(0)
        record = f.read(1024)

        self.locfmt = record[88:96]

        if self.locfmt == b'BIG-IEEE':
            self.endian = '>'
        elif self.locfmt == b'LTL-IEEE':
            self.endian = '<'
        else:
            raise ValueError('unrecognized format: {0!r}'.format(self.locfmt))

        (self.locidw, self.nd, self.ni, self.locifn, self.fward, self.bward,
         self.free) = struct.unpack(self.endian + '8sII60sIII', record[:88])

        if self.locidw != b'DAF/SPK ':
            raise ValueError(
                'the first bytes do not identify this as a DAF/SPK file')

        print vars(self)


def main():
    with open('jup310.tmp', 'rb') as f:
        daf = DAF(f)

    # if b[500:1000].strip(b'\0') != FTPSTR:
    #     raise ValueError('the file has been damaged')

    # print repr(b[128:128+32])
    # print b.index('LTL-IEEE')

main()
