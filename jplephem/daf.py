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
        if getattr(file_object, 'encoding', None):
            raise ValueError('file_object must be opened in binary "b" mode')

        self.file = file_object
        self.fileno = file_object.fileno()

        file_record = self.read_record(1)

        self.locidw = file_record[:8].upper().rstrip()
        if not self.locidw.startswith(b'DAF/'):
            raise ValueError('file starts with {0!r}, not the 4 bytes {0!r}'
                             .format(self.locidw, b'DAF/'))

        if file_record[500:1000].strip(b'\0') != FTPSTR:
            raise ValueError('this SPK file has been damaged')

        self.locfmt = file_record[88:96]
        self.endian = LOCFMT.get(self.locfmt)
        if self.endian is None:
            raise ValueError('unsupported format {0!r}'.format(self.locfmt))

        (self.nd, self.ni, locifn, self.fward, self.bward, self.free
         ) = struct.unpack(self.endian + 'II60sIII', file_record[8:88])

        summary_size = self.nd + (self.ni + 1) // 2
        self.summary_step = 8 * summary_size
        self.summary_format = self.endian + 'd' * self.nd + 'i' * self.ni
        self.summary_length = struct.calcsize(self.summary_format)
        self.locifn = locifn.upper().rstrip()

    def read_record(self, n):
        """Return record `n` as 1,024 bytes; records are indexed from 1."""
        self.file.seek(n * K - 1024)
        return self.file.read(K)

    def map_words(self, start, end):
        """Return a memory-map of the elements `start` through `end`.

        The memory map will offer the 8-byte double-precision floats
        ("elements") in the file from index `start` through to the index
        `end`, inclusive, both counting the first float as element 1.
        Memory maps must begin on a page boundary, so `skip` returns the
        number of extra bytes at the beginning of the return value.

        """
        i, j = 8 * start - 8, 8 * end
        skip = i % mmap.PAGESIZE
        r = mmap.ACCESS_READ
        m = mmap.mmap(self.fileno, length=j-i+skip, access=r, offset=i-skip)
        if sys.version_info > (3,):
            m = memoryview(m)  # so further slicing can return views
        return m, skip

    def comments(self):
        """Return the text inside the comment area of the file."""
        record_numbers = range(2, self.fward)
        if not record_numbers:
            return ''
        data = b''.join(self.read_record(n)[0:1000] for n in record_numbers)
        try:
            return data[:data.find(b'\4')].decode('ascii').replace('\0', '\n')
        except IndexError:
            raise ValueError('DAF file comment area is missing its EOT byte')
        except UnicodeDecodeError:
            raise ValueError('DAF file comment area is not ASCII text')

    def map_array(self, start, end):
        """Return floats from `start` to `end` inclusive, indexed from 1."""
        data, skip = self.map_words(start, end)
        skip //= 8
        return ndarray(end - start + 1 + skip, self.endian + 'd', data)[skip:]

    def summaries(self):
        """Yield (name, (value, value, ...)) for each summary in the file."""

        record_number = self.fward
        length = self.summary_length
        step = self.summary_step

        while record_number:
            summary_record = self.read_record(record_number)
            name_record = self.read_record(record_number + 1)

            next_number, previous_number, n_summaries = struct.unpack(
                self.endian + 'ddd', summary_record[:24])

            for i in range(0, int(n_summaries) * step, step):
                j = i + 24
                name = name_record[i:i+step].strip()
                data = summary_record[j:j+length]
                values = struct.unpack(self.summary_format, data)
                yield name, values

            record_number = int(next_number)
