# taken from research/kzg_data_availability

from py_ecc import optimized_bls12_381 as b


def _simple_ft(vals, modulus, roots_of_unity):
    L = len(roots_of_unity)
    o = []
    for i in range(L):
        last = b.Z1 if type(vals[0]) == tuple else 0
        for j in range(L):
            if type(vals[0]) == tuple:
                last = b.add(last, b.multiply(vals[j], roots_of_unity[(i*j)%L]))
            else:
                last += vals[j] * roots_of_unity[(i*j)%L]
        o.append(last if type(last) == tuple else last % modulus)
    return o

def _fft(vals, modulus, roots_of_unity):
    if len(vals) <= 4 and type(vals[0]) != tuple:
        #return vals
        return _simple_ft(vals, modulus, roots_of_unity)
    elif len(vals) == 1 and type(vals[0]) == tuple:
        return vals
    L = _fft(vals[::2], modulus, roots_of_unity[::2])
    R = _fft(vals[1::2], modulus, roots_of_unity[::2])
    o = [0 for i in vals]
    for i, (x, y) in enumerate(zip(L, R)):
        y_times_root = b.multiply(y, roots_of_unity[i]) if type(y) == tuple else y*roots_of_unity[i]
        o[i] = b.add(x, y_times_root) if type(x) == tuple else (x+y_times_root) % modulus
        o[i+len(L)] = b.add(x, b.neg(y_times_root)) if type(x) == tuple else (x-y_times_root) % modulus
    return o


def expand_root_of_unity(root_of_unity, modulus):
    # Build up roots of unity
    rootz = [1, root_of_unity]
    while rootz[-1] != 1:
        rootz.append((rootz[-1] * root_of_unity) % modulus)
    return rootz


def fft(vals, modulus, root_of_unity, inv=False):
    rootz = expand_root_of_unity(root_of_unity, modulus)
    # Fill in vals with zeroes if needed
    if len(rootz) > len(vals) + 1:
        vals = vals + [0] * (len(rootz) - len(vals) - 1)
    if inv:
        # Inverse FFT
        invlen = pow(len(vals), modulus-2, modulus)
        if type(vals[0]) == tuple:
            return [b.multiply(x, invlen) for x in
                    _fft(vals, modulus, rootz[:0:-1])]
        else:
            return [(x*invlen) % modulus for x in
                    _fft(vals, modulus, rootz[:0:-1])]
    else:
        # Regular FFT
        return _fft(vals, modulus, rootz[:-1])
