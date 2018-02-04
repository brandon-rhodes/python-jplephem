"""Access a NASA JPL SPICE Double Precision Array File (DAF).

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html

"""
import mmap
import sys
from struct import Struct
from numpy import ndarray

FTPSTR = b'FTPSTR:\r:\n:\r\n:\r\x00:\x81:\x10\xce:ENDFTP'  # FTP test string
LOCFMT = {b'BIG-IEEE': '>', b'LTL-IEEE': '<'}
K = 1024

class DAF(object):
    """Access to NASA SPICE Double Precision Array Files (DAF).

    Provide the constructor with a ``file_object`` for full access to
    both the segment summaries and to the numeric arrays.  If you pass a
    ``StringIO`` instead, then you can fetch the summary information but
    not access the arrays.

    """
    def __init__(self, file_object):
        if getattr(file_object, 'encoding', None):
            raise ValueError('file_object must be opened in binary "b" mode')

        self.file = file_object

        file_record = self.read_record(1)

        def unpack():
            fmt = self.endian + '8sII60sIII8s603s28s297s'
            self.file_record_struct = Struct(fmt)
            (locidw, self.nd, self.ni, self.locifn, self.fward, self.bward,
             self.free, locfmt, self.prenul, self.ftpstr, self.pstnul
            ) = self.file_record_struct.unpack(file_record)

        self.locidw = file_record[:8].upper().rstrip()

        if self.locidw == b'NAIF/DAF':
            for self.locfmt, self.endian in LOCFMT.items():
                unpack()
                if self.nd == 2:
                    break
            else:
                raise ValueError('neither a big- nor a little-endian scan'
                                 ' of this file produces the expected ND=2')
        elif self.locidw.startswith(b'DAF/'):
            if file_record[500:1000].strip(b'\0') != FTPSTR:
                raise ValueError('this SPK file has been damaged')
            self.locfmt = file_record[88:96]
            self.endian = LOCFMT.get(self.locfmt)
            if self.endian is None:
                raise ValueError('unknown format {0!r}'.format(self.locfmt))
            unpack()
        else:
            raise ValueError('file starts with {0!r}, not "NAIF/DAF" or "DAF/"'
                             .format(self.locidw))

        fmt = self.endian + 'd' * self.nd + 'i' * self.ni

        self.summary_struct = struct = Struct(fmt)
        self.summary_length = length = struct.size
        self.summary_step = length + (-length % 8) # pad to 8 bytes
        self.locifn_text = self.locifn.rstrip()

    def read_record(self, n):
        """Return record `n` as 1,024 bytes; records are indexed from 1."""
        self.file.seek(n * K - 1024)
        return self.file.read(K)

    def write_file_record(self):
        data = self.file_record_struct.pack(
            self.locidw, self.nd, self.ni, self.locifn, self.fward, self.bward,
            self.free, self.locfmt, self.prenul, self.ftpstr, self.pstnul,
        )
        self.file.seek(0)
        self.file.write(data)

    def map_words(self, start, end):
        """Return a memory-map of the elements `start` through `end`.

        The memory map will offer the 8-byte double-precision floats
        ("elements") in the file from index `start` through to the index
        `end`, inclusive, both counting the first float as element 1.
        Memory maps must begin on a page boundary, so `skip` returns the
        number of extra bytes at the beginning of the return value.

        """
        fileno = self.file.fileno() # requires a true file object
        i, j = 8 * start - 8, 8 * end
        skip = i % mmap.ALLOCATIONGRANULARITY
        r = mmap.ACCESS_READ
        m = mmap.mmap(fileno, length=j-i+skip, access=r, offset=i-skip)
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

    def read_array(self, start, end):
        """Return floats from `start` to `end` inclusive, indexed from 1.

        The entire range of floats is immediately read into memory from
        the file, making this efficient for small sequences of floats
        whose values are all needed immediately.

        """
        f = self.file
        f.seek(8 * (start - 1))
        length = 1 + end - start
        data = f.read(8 * length)
        return ndarray(length, self.endian + 'd', data)

    def map_array(self, start, end):
        """Return floats from `start` to `end` inclusive, indexed from 1.

        Instead of pausing to load all of the floats into RAM, this
        routine creates a memory map which will load data from the file
        only as it is accessed, and then will let it expire back out to
        disk later.  This is very efficient for large data sets to which
        you need random access.

        """
        data, skip = self.map_words(start, end)
        skip //= 8
        return ndarray(end - start + 1 + skip, self.endian + 'd', data)[skip:]

    def summary_records(self):
        """Yield (record_number, n_summaries, summary_data) for each record.

        Readers will only use the second two values in each tuple.
        Writers can update the record using the `record_number`.

        """
        record_number = self.fward
        unpack = Struct(self.endian + 'ddd').unpack
        while record_number:
            data = self.read_record(record_number)
            next_number, previous_number, n_summaries = unpack(data[:24])
            yield record_number, n_summaries, data
            record_number = int(next_number)

    def summaries(self):
        """Yield (name, (value, value, ...)) for each summary in the file."""
        length = self.summary_length
        step = self.summary_step
        for record_number, n_summaries, summary_data in self.summary_records():
            name_data = self.read_record(record_number + 1)
            for i in range(0, int(n_summaries) * step, step):
                j = 8*3 + i
                name = name_data[i:i+step].strip()
                data = summary_data[j:j+length]
                values = self.summary_struct.unpack(data)
                yield name, values


NAIF_DAF = DAF  # a separate class supported NAIF/DAF format in jplephem 2.2
