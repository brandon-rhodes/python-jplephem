"""Access a NASA JPL SPICE Double Precision Array File (DAF).

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html

"""
import io
import mmap
import sys
try:  # Use low-level module to avoid importing huge 'threading.py'
    from _thread import allocate_lock
except:
    from thread import allocate_lock
from struct import Struct
from numpy import array as numpy_array, ndarray

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
        self.lock = allocate_lock()
        self._map = None
        self._array = None

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

        self.locifn_text = self.locifn.rstrip()

        summary_format = 'd' * self.nd + 'i' * self.ni

        self.summary_control_struct = Struct(self.endian + 'ddd')
        self.summary_struct = struct = Struct(self.endian + summary_format)
        self.summary_length = length = struct.size
        self.summary_step = length + (-length % 8) # pad to 8 bytes
        self.summaries_per_record = (1024 - 8 * 3) // self.summary_step

    def read_record(self, n):
        """Return record `n` as 1,024 bytes; records are indexed from 1."""
        with self.lock:
            self.file.seek(n * K - K)
            return self.file.read(K)

    def write_record(self, n, data):
        """Write `data` to file record `n`; records are indexed from 1."""
        with self.lock:
            self.file.seek(n * K - K)
            return self.file.write(data)

    def write_file_record(self):
        data = self.file_record_struct.pack(
            self.locidw.ljust(8, b' '), self.nd, self.ni, self.locifn,
            self.fward, self.bward, self.free, self.locfmt,
            self.prenul, self.ftpstr, self.pstnul,
        )
        self.write_record(1, data)

    def map_words(self, start, end):
        """Return a memory-map of the elements `start` through `end`.

        The memory map will offer the 8-byte double-precision floats
        ("elements") in the file from index `start` through to the index
        `end`, inclusive, both counting the first float as element 1.
        Memory maps must begin on a page boundary, so `skip` returns the
        number of extra bytes at the beginning of the return value.

        If a memory map is not available on your operating system, then
        the segment's bytes are simply read into an array instead.

        """
        i, j = 8 * start - 8, 8 * end
        try:
            fileno = self.file.fileno()
        except (AttributeError, io.UnsupportedOperation):
            m = None
        else:
            skip = i % mmap.ALLOCATIONGRANULARITY
            r = mmap.ACCESS_READ
            try:
                m = mmap.mmap(fileno, length=j-i+skip, access=r, offset=i-skip)
            except OSError:
                m = None
        if m is None:
            skip = 0
            with self.lock:
                self.file.seek(i)
                m = self.file.read(j - i)
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
        length = 1 + end - start
        with self.lock:
            f.seek(8 * (start - 1))
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
        if self._array is None:
            self._map, skip = self.map_words(1, self.free - 1)
            assert skip == 0
            self._array = ndarray(self.free - 1, self.endian + 'd', self._map)
        return self._array[start - 1 : end]

    def summary_records(self):
        """Yield (record_number, n_summaries, record_data) for each record.

        Readers will only use the second two values in each tuple.
        Writers can update the record using the `record_number`.

        """
        record_number = self.fward
        unpack = self.summary_control_struct.unpack
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
                j = self.summary_control_struct.size + i
                name = name_data[i:i+step].strip()
                data = summary_data[j:j+length]
                values = self.summary_struct.unpack(data)
                yield name, values

    def map(self, summary_values):
        """Return the array of floats described by a summary.

        Instead of pausing to load all of the floats into RAM, this
        routine creates a memory map which will load data from the file
        only as it is accessed, and then will let it expire back out to
        disk later.  This is very efficient for large data sets to which
        you need random access.

        """
        return self.map_array(summary_values[-2], summary_values[-1])

    def add_array(self, name, values, array):
        """Add a new array to the DAF file.

        The summary will be initialized with the `name` and `values`,
        and will have its start word and end word fields set to point to
        where the `array` of floats has been appended to the file.

        This method is not thread-safe.

        """
        f = self.file
        scs = self.summary_control_struct

        record_number = self.bward
        data = bytearray(self.read_record(record_number))
        next_record, previous_record, n_summaries = scs.unpack(data[:24])

        if n_summaries < self.summaries_per_record:
            summary_record = record_number
            name_record = summary_record + 1
            data[:24] = scs.pack(next_record, previous_record, n_summaries + 1)
            self.write_record(summary_record, data)
        else:
            summary_record = ((self.free - 1) * 8 + 1023) // 1024 + 1
            name_record = summary_record + 1
            free_record = summary_record + 2

            data[:24] = scs.pack(summary_record, previous_record, n_summaries)
            self.write_record(record_number, data)

            n_summaries = 0
            summaries = scs.pack(0, record_number, 1).ljust(1024, b'\0')
            names = b'\0' * 1024
            self.write_record(summary_record, summaries)
            self.write_record(name_record, names)

            self.bward = summary_record
            self.free = (free_record - 1) * 1024 // 8 + 1

        array = numpy_array(array, self.endian + 'f8')

        start_word = self.free
        f.seek((start_word - 1) * 8)
        f.write(array.view())
        end_word = f.tell() // 8

        self.free = end_word + 1
        self.write_file_record()

        values = values[:self.nd + self.ni - 2] + (start_word, end_word)

        base = 1024 * (summary_record - 1)
        offset = int(n_summaries) * self.summary_step
        f.seek(base + scs.size + offset)
        f.write(self.summary_struct.pack(*values))
        f.seek(base + 1024 + offset)
        f.write(name[:self.summary_length].ljust(self.summary_step, b' '))


NAIF_DAF = DAF  # a separate class supported NAIF/DAF format in jplephem 2.2
