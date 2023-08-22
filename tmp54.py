import numpy as np
from sys import argv
from jplephem.names import target_names
from jplephem.spk import SPK, titlecase

_component_names = {1: 'x', 2: 'x,y', 3: 'x,y,z'}

def _compute_calendar_date(jd_integer, julian_before=None):
    """Convert Julian day ``jd_integer`` into a calendar (year, month, day).

    Uses the proleptic Gregorian calendar unless ``julian_before`` is
    set to a specific Julian day, in which case the Julian calendar is
    used for dates older than that.

    """
    use_gregorian = (julian_before is None) or (jd_integer >= julian_before)

    # See the Explanatory Supplement to the Astronomical Almanac 15.11.
    f = jd_integer + 1401
    f += use_gregorian * ((4 * jd_integer + 274277) // 146097 * 3 // 4 - 38)
    e = 4 * f + 3
    g = e % 1461 // 4
    h = 5 * g + 2
    day = h % 153 // 5 + 1
    month = (h // 153 + 2) % 12 + 1
    year = e // 1461 - 4716 + (12 + 2 - month) // 12
    return year, month, day

def _format_date(jd):
    year, month, day = _compute_calendar_date(int(jd))
    return '{:4}-{:02}-{:02}'.format(year, month, day)

def print_type2_segment(segment):
    initial_epoch, interval_length, coefficients = segment.load_array()
    component_count, n, coefficient_count = coefficients.shape

    center = titlecase(target_names.get(segment.center, 'Unknown center'))
    target = titlecase(target_names.get(segment.target, 'Unknown target'))

    start_date = _format_date(segment.start_jd)
    end_date = _format_date(segment.end_jd)
    print(f'Date: {start_date} to {end_date}')
    print(f'Center: {center} ({segment.center})')
    print(f'Target: {target} ({segment.target})')

    c = _component_names.get(component_count, str(component_count))
    s = 8 * (segment.end_i - segment.start_i)
    print(f'Components: {c}  '
          f'Coefficient count: {coefficient_count}  '
          f'Size: {s:,} bytes')

    indent = '\n' + ' ' * 14
    if n > 5:
        ii = [0, 1, 2, '   ...', n-3, n-2, n-1]
    else:
        ii = range(n)

    with np.printoptions(precision=1, linewidth=82 - len(indent)):
        for i in ii:
            if isinstance(i, str):
                print(i)
                continue
            d = initial_epoch + i * interval_length
            nums = str(coefficients[:,i]).replace('\n', indent)
            print(_format_date(d).rjust(13), nums)

        i += 1
        d = initial_epoch + i * interval_length
        print(_format_date(d).rjust(13))

if __name__ == '__main__':
    kernel = SPK.open(argv[1])  # 'de441.bsp' or 'de441_partial.bsp'
    segment_index = int(argv[2])  # like '2' or '0'
    print_type2_segment(kernel.segments[segment_index])
