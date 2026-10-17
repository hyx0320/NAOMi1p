import struct, zlib, sys, os
import numpy as np

# ---------- minimal TIFF reader: classic + BigTIFF, strips, none/packbits/deflate ----------
TSIZE = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2, 9: 4, 10: 8, 11: 4, 12: 8}
TCODE = {1: 'B', 3: 'H', 4: 'I', 6: 'b', 8: 'h', 9: 'i', 11: 'f', 12: 'd'}

def _read_ifds(path):
    with open(path, 'rb') as f:
        head = f.read(16)
    en = '<' if head[:2] == b'II' else '>'
    magic = struct.unpack(en + 'H', head[2:4])[0]
    if magic == 42:
        big, off = False, struct.unpack(en + 'I', head[4:8])[0]
    elif magic == 43:
        big, off = True, struct.unpack(en + 'Q', head[8:16])[0]
    else:
        raise ValueError('not a TIFF')
    ifds = []
    with open(path, 'rb') as f:
        while off:
            f.seek(off)
            if big:
                n = struct.unpack(en + 'Q', f.read(8))[0]
                ents = [struct.unpack(en + 'HHQQ', f.read(20)) for _ in range(n)]
                nxt = struct.unpack(en + 'Q', f.read(8))[0]
            else:
                n = struct.unpack(en + 'H', f.read(2))[0]
                ents = [struct.unpack(en + 'HHII', f.read(12)) for _ in range(n)]
                nxt = struct.unpack(en + 'I', f.read(4))[0]
            tags = {}
            for tag, typ, cnt, val in ents:
                sz = TSIZE.get(typ, 1) * cnt
                inline = 8 if big else 4
                if sz <= inline:
                    codes = TCODE.get(typ)
                    if codes is None:
                        vals = [val]
                    else:
                        nb = TSIZE[typ] * cnt
                        raw = val.to_bytes(inline, 'little' if en == '<' else 'big')[:nb]
                        vals = list(struct.unpack(en + str(cnt) + codes, raw))
                else:
                    pos = f.tell()
                    f.seek(val)
                    codes = TCODE.get(typ)
                    raw = f.read(sz)
                    vals = list(struct.unpack(en + str(cnt) + codes, raw)) if codes else []
                    f.seek(pos)
                tags[tag] = vals
            ifds.append(tags)
            off = nxt
    return ifds, en

def _unpackbits(raw, n):
    out = bytearray(); i = 0
    while i < len(raw) and len(out) < n:
        c = raw[i]; i += 1
        if c < 128:
            out += raw[i:i + c + 1]; i += c + 1
        elif c > 128:
            if i < len(raw):
                out += bytes([raw[i]]) * (257 - c)
            i += 1
    return bytes(out)

def read_tiff(path, max_pages=None):
    ifds, en = _read_ifds(path)
    pages = []
    with open(path, 'rb') as f:
        for t in ifds[:max_pages] if max_pages else ifds:
            w, h = t[256][0], t[257][0]
            bps, comp = t.get(258, [8])[0], t.get(259, [1])[0]
            spp = t.get(277, [1])[0]
            offs = t[273]; cnts = t[279]
            nbytes = w * h * (bps // 8) * spp
            buf = bytearray()
            for o, c in zip(offs, cnts):
                f.seek(o)
                raw = f.read(c)
                if comp == 1:
                    buf += raw
                elif comp == 32773:
                    buf += _unpackbits(raw, w * h * (bps // 8) * spp)
                elif comp == 8:
                    buf += zlib.decompress(raw)
                else:
                    raise ValueError('unsupported compression %s' % comp)
            buf = bytes(buf[:nbytes]).ljust(nbytes, b'\x00')
            dt = en + ('u1' if bps == 8 else 'u2' if bps == 16 else 'u4')
            arr = np.frombuffer(buf, dtype=dt).reshape(h, w).astype(np.float32)
            pages.append(arr)
    return np.stack(pages)

def norm(a, lo=1.0, hi=99.5):
    v = a[np.isfinite(a)]
    if v.size == 0:
        return np.zeros_like(a)
    p1, p2 = np.percentile(v, [lo, hi])
    return np.clip((a - p1) / max(p2 - p1, 1e-9), 0, 1)

def write_png(path, arr):
    h, w = arr.shape
    raw = b''.join(b'\x00' + arr[y].tobytes() for y in range(h))
    def chunk(tag, data):
        return (struct.pack('>I', len(data)) + tag + data +
                struct.pack('>I', zlib.crc32(tag + data) & 0xffffffff))
    png = (b'\x89PNG\r\n\x1a\n'
           + chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 0, 0, 0, 0))
           + chunk(b'IDAT', zlib.compress(raw, 6)) + chunk(b'IEND', b''))
    with open(path, 'wb') as f:
        f.write(png)

R = r'C:\Users\YaXinHu\Desktop\NAOMi1p'
O = r'D:\MesoTwins\NAOMi1p_output\w_dend\res_0.80\vol_600_150_NA_0.35_Hz_10_exp_5_d_30k_pw_0.50\11'

print('neur_vol ...', end=' ', flush=True); vol = read_tiff(os.path.join(O, 'neur_vol.tiff'))
print(vol.shape)
print('vessels  ...', end=' ', flush=True); ves = read_tiff(os.path.join(O, 'neur_ves_all.tiff'))
print(ves.shape)
print('psf      ...', end=' ', flush=True); psf = read_tiff(os.path.join(O, 'psf.tiff'))
print(psf.shape)
print('mask     ...', end=' ', flush=True); msk = read_tiff(os.path.join(O, 'psf_mask.tiff'))
print(msk.shape)
print('frame    ...', end=' ', flush=True); frm = read_tiff(os.path.join(R, r'TimeTraceCode\sim_frame_0001.tif'))
print(frm.shape)

S = 300
tiles = [
    ('volume z=30', vol[30][::2, ::2]),
    ('volume z=75', vol[75][::2, ::2]),
    ('volume z=120', vol[120][::2, ::2]),
    ('volume z-MIP', vol.max(0)[::2, ::2]),
    ('vasculature z=100', ves[100][::3, ::3][:S, :S]),
    ('PSF mask / illumination', msk[0][::2, ::2]),
    ('PSF xz projection', np.repeat(np.repeat(psf.max(1), 8, 0), 8, 1)),
    ('exported slice sim_frame_0001', frm[0][::2, ::2]),
]
canvas = np.zeros((2 * S, 4 * S), np.uint8)
for i, (name, t) in enumerate(tiles):
    r, c = divmod(i, 4)
    t = norm(t)
    h, w = min(t.shape[0], S), min(t.shape[1], S)
    canvas[r*S:r*S+h, c*S:c*S+w] = (t[:h, :w] * 255).astype(np.uint8)

out = os.path.join(R, r'.mlenv\preview_volume.png')
write_png(out, canvas)
print('wrote', out, canvas.shape)
print('volume: min=%.0f max=%.0f, nonzero=%.2f%%' % (vol.min(), vol.max(), 100.0*(vol > 0).mean()))
print('vessels: max=%.0f  psf: max=%.3f' % (ves.max(), psf.max()))
