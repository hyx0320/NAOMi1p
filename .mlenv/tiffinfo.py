import struct, os, sys

def pages_info(p):
    with open(p, 'rb') as f:
        b = f.read(16)
    if b[:2] not in (b'II', b'MM'):
        return 'not a TIFF: %r' % b[:4]
    en = '<' if b[:2] == b'II' else '>'
    magic = struct.unpack(en + 'H', b[2:4])[0]
    if magic == 42:
        big = False
        off = struct.unpack(en + 'I', b[4:8])[0]
    elif magic == 43:
        big = True
        off = struct.unpack(en + 'Q', b[8:16])[0]
    else:
        return 'unknown TIFF magic %d' % magic

    bpsmap = {1: 'uint8', 2: 'uint16', 3: 'int16', 4: 'uint32', 5: 'int32',
              6: 'float32', 7: 'float64', 11: 'float32'}
    compmap = {1: 'none', 5: 'lzw', 7: 'jpeg', 8: 'deflate', 32773: 'packbits'}

    with open(p, 'rb') as f:
        ifd = off
        page = 0
        first = None
        while ifd:
            f.seek(ifd)
            if big:
                n = struct.unpack(en + 'Q', f.read(8))[0]
                ents = [struct.unpack(en + 'HHQQ', f.read(20)) for _ in range(n)]
                ifd = struct.unpack(en + 'Q', f.read(8))[0]
            else:
                n = struct.unpack(en + 'H', f.read(2))[0]
                ents = [struct.unpack(en + 'HHII', f.read(12)) for _ in range(n)]
                ifd = struct.unpack(en + 'I', f.read(4))[0]
            tags = {}
            for tag, typ, cnt, val in ents:
                tags[tag] = val
            if page == 0:
                first = (tags.get(256), tags.get(257), tags.get(258), tags.get(259))
            page += 1
            if page > 5000:
                break
    w, h, bps, comp = first
    return '%sx%s  %s  compression=%s  pages=%d' % (
        w, h, bpsmap.get(bps, 'bps=%s' % bps), compmap.get(comp, comp), page)

for p in sys.argv[1:]:
    if not os.path.exists(p):
        print('MISSING', p)
        continue
    try:
        print('%9.2f MB  %-28s %s' % (os.path.getsize(p) / 1048576.0,
                                      os.path.basename(p), pages_info(p)))
    except Exception as e:
        print('ERROR', p, e)
