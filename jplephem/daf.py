"""Access a NASA JPL SPICE Double Precision Array File (DAF).

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html

"""
import mmap
import struct
import sys
from numpy import ndarray

FTPSTR = b'FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP'  # FTP test string
LOCFMT = {b'BIG-IEEE': '>', b'LTL-IEEE': '<'}
K = 1024

class DAF(object):
    """Access to NASA SPICE Double Precision Array Files (DAF)."""

    def __init__(self, file_object):
        self.map = mmap.mmap(file_object.fileno(), 0, access=mmap.ACCESS_READ)
        if sys.version_info > (3,):
            self.map = memoryview(self.map)

        self.locidw = bytes(self.map[:8]).upper().rstrip()
        if not self.locidw.startswith(b'DAF/'):
            raise ValueError('file starts with {0!r}, not the 4 bytes {0!r}'
                             .format(self.locidw, b'DAF/'))

        if bytes(self.map[500:1000]).strip(b'\0') != FTPSTR:
            raise ValueError('this SPK file has been damaged')

        self.locfmt = bytes(self.map[88:96])
        self.endian = LOCFMT.get(self.locfmt)
        if self.endian is None:
            raise ValueError('unsupported format {0!r}'.format(self.locfmt))

        (self.nd, self.ni, locifn, self.fward, self.bward, self.free
         ) = struct.unpack(self.endian + 'II60sIII', self.map[8:88])

        summary_size = self.nd + (self.ni + 1) // 2
        self.summary_step = 8 * summary_size
        self.summary_format = self.endian + 'd' * self.nd + 'i' * self.ni
        self.summary_length = struct.calcsize(self.summary_format)
        self.locifn = locifn.upper().rstrip()

    def bytes(self, start, end):
        """Return data of words `start` to `end` inclusive; indexed from 1."""
        return self.map[8 * start - 8 : 8 * end]

    def record(self, n):
        """Return record `n` as bytes; records are indexed from 1."""
        stop = n * K
        return self.map[stop - K:stop]

    def comments(self):
        """Return the text inside the comment area of the file."""
        record_numbers = range(2, self.fward)
        if not record_numbers:
            return ''
        data = b''.join(bytes(self.record(n)[0:1000]) for n in record_numbers)
        try:
            return data[:data.find(b'\4')].decode('ascii').replace('\0', '\n')
        except IndexError:
            raise ValueError('DAF file comment area is missing its EOT byte')
        except UnicodeDecodeError:
            raise ValueError('DAF file comment area is not ASCII text')

    def array(self, start, end):
        """Return floats from `start` to `end` inclusive; indexed from 1."""
        data = self.bytes(start, end)
        return ndarray(end - start + 1, self.endian + 'd', data)

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
