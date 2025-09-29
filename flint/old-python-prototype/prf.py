# PRF that outputs a number in [0, t)
import hashlib


def prf(t: int, key: bytes, data: bytes) -> int:
    # Use + 16 to ensure the output distribution is statistically close to uniform over [0, t]
    # 16 * 8 = 128, so 2^{-128} is the statistical distance
    # TODO: can probably make this faster with AES vs SHAKE
    return int.from_bytes(hashlib.shake_256(key + data).digest(int(t).bit_length() // 8 + 16), byteorder='big') % t

if __name__ == "__main__":
    # Test PRF
    for t in [10, 100, 10000, 2 ** 256, 2 ** 512]:
        print(f't =             {t}')
        for data in [b'data0', b'data1', b'data2', b'data3']:
            print(f'PRF({data}) =', prf(t, ('key' + str(t)).encode(), data))
    assert prf(100, b'key', b'data') == prf(100, b'key', b'data')
