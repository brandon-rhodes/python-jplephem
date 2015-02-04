"""Access a NASA JPL SPICE Double Precision Array File (DAF).

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html

"""
import mmap
import struct
import sys

BFF = b'BIG-IEEE', b'LTL-IEEE', b'VAX-GFLT', b'VAX-DFLT'  # Binary file formats
FTPSTR = b'FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP' # FTP test string
RECORD_LENGTH = 1024

class DAF(object):
    """Access to NASA SPICE Double Precision Array Files (DAF)."""

    def __init__(self, path):
        with open(path, 'rb') as f:
            m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        if sys.version_info > (3,):
            m = memoryview(m)
        self.map = m

        self.locfmt = m[88:96]
        if self.locfmt == b'BIG-IEEE':
            self.endian = '>'
        elif self.locfmt == b'LTL-IEEE':
            self.endian = '<'
        else:
            raise ValueError('unsupported format {0!r}'.format(self.locfmt))

        (locidw, nd, ni, locifn, self.fward, self.bward, self.free
         ) = struct.unpack(self.endian + '8sII60sIII', m[:88])

        summary_size = nd + (ni + 1) // 2
        self.summary_step = 8 * summary_size
        self.summary_format = self.endian + 'd' * nd + 'i' * ni
        self.summary_length = struct.calcsize(self.summary_format)

        self.locidw = locidw.upper().rstrip()
        self.locifn = locifn.upper().rstrip()
        if self.locidw != b'DAF/SPK':
            raise ValueError('the first file bytes do not specify "DAF/SPK"')

        if bytes(m[500:1000]).strip(b'\0') != FTPSTR:
            raise ValueError('this SPK file has been damaged')

    def bytes(self, start, stop):
        """Return data from word `start` to `stop`, inclusive."""
        return self.map[8 * start - 8 : 8 * stop]

    def record(self, n):
        """Return record `n` as bytes; records are indexed from 1."""
        stop = RECORD_LENGTH * n
        return self.map[stop - RECORD_LENGTH:stop]

    def summaries(self):
        """Yield (name, (value, value, ...)) for each summary in the file."""

        record_number = self.fward
        length = self.summary_length
        step = self.summary_step

        while record_number:
            summary_record = self.record(record_number)
            name_record = self.record(record_number + 1)

            next_number, previous_number, n_summaries = struct.unpack(
                self.endian + 'ddd', summary_record[:24])

            for i in range(0, int(n_summaries) * step, step):
                j = i + 24
                name = bytes(name_record[i:i+step]).strip()
                data = summary_record[j:j+length]
                values = struct.unpack(self.summary_format, data)
                yield name, values

            record_number = int(next_number)
