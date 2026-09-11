import struct

import pytest

from gsynth_engine import chromatogram
from gsynth_engine.sequence import SequenceError
from gsynth_engine.tests.test_chromatogram import READ, build_ab1, build_scf


@pytest.mark.parametrize('reader,data,message', [
    (chromatogram.read_scf, b'.scf', 'too small'),
    (chromatogram.read_scf, bytes(128), 'not an SCF'),
    (chromatogram.read_trace, b'junk', 'too small'),
])
def test_invalid_trace_headers(reader, data, message):
    with pytest.raises(SequenceError, match=message):
        reader(data)


@pytest.mark.parametrize('offset,value,message', [
    (36, b'bad!', 'invalid format version'),
    (36, b'\xff.00', 'invalid format version'),
    (36, b'4.00', 'not supported'),
    (4, struct.pack('>I', 0), 'samples'),
    (4, struct.pack('>I', chromatogram.MAX_SAMPLES + 1), 'samples'),
    (12, struct.pack('>I', 0), 'base calls'),
    (12, struct.pack('>I', chromatogram.MAX_BASES + 1), 'base calls'),
])
def test_scf_rejects_invalid_version_and_counts(offset, value, message):
    data = bytearray(build_scf(READ))
    data[offset:offset + 4] = value
    with pytest.raises(SequenceError, match=message):
        chromatogram.read_scf(bytes(data))


def test_scf_private_section_is_bounds_checked():
    data = bytearray(build_scf(READ))
    data[48:52] = struct.pack('>I', 4)
    with pytest.raises(SequenceError, match='private data'):
        chromatogram.read_scf(bytes(data))
    trace = chromatogram.read_scf(bytes(data) + b'meta')
    assert trace.sequence == READ


def test_trace_empty_window_and_single_peak():
    trace = chromatogram.read_ab1(build_ab1('A'))
    assert trace._spacing() == 12
    window = trace.alignment_track(0, 0)
    assert window['sequence'] == ''
    assert window['sample_count'] == 0
    assert all(values == [] for values in window['traces'].values())


@pytest.mark.parametrize('kind,payload,expected', [
    (18, b'\x03abc', 'abc'), (19, b'abc\x00', 'abc'), (99, b'abc\x00', b'abc\x00'),
])
def test_abif_inline_string_and_opaque_values(kind, payload, expected):
    entry = chromatogram._read_entry(struct.pack('>4sihhii4si', b'TEST', 1, kind, 1, 4, 4, payload, 0), 0)
    assert chromatogram._value(b'', entry) == expected


@pytest.mark.parametrize('kind', [2, 99])
def test_abif_quality_byte_encodings(kind):
    data = bytearray(build_ab1(READ))
    data[128 + 28 + 8:128 + 28 + 10] = struct.pack('>h', kind)
    assert chromatogram.read_ab1(bytes(data)).quality == [40] * len(READ)


def test_abif_truncated_directory_is_rejected():
    data = bytearray(build_ab1(READ))
    data[18:22] = struct.pack('>i', 10000)
    with pytest.raises(SequenceError):
        chromatogram.read_ab1(bytes(data[:140]))


def test_abif_signal_limits(monkeypatch):
    data = build_ab1(READ)
    monkeypatch.setattr(chromatogram, 'MAX_SAMPLES', 1)
    with pytest.raises(SequenceError, match='samples'):
        chromatogram.read_ab1(data)


def test_abif_peaks_beyond_signal_are_rejected():
    data = build_ab1(READ, peaks=[10000] * len(READ))
    with pytest.raises(SequenceError, match='stops before its last base'):
        chromatogram.read_ab1(data)
