class BitfieldRow:
    def __init__(self, row_tuple):
        # Accepts a tuple/list of ints and converts to bitfield
        bitfield = 0
        for k in row_tuple:
            bitfield |= 1 << (k - 1)
        self.bitfield = bitfield

    def __contains__(self, item):
        return (self.bitfield & (1 << (item - 1))) != 0

    def __iter__(self):
        # Iterate over set bits in descending order
        bits = []
        bf = self.bitfield
        k = 1
        while bf:
            if bf & 1:
                bits.append(k)
            bf >>= 1
            k += 1
        yield from reversed(bits)

    def __len__(self):
        return bin(self.bitfield).count("1")

    def __getitem__(self, idx):
        bits = [i for i in self]  # noqa: C416
        return bits[idx]

    def to_tuple(self):
        return tuple(self)

    def __repr__(self):
        return f"BitfieldRow({self.to_tuple()})"

    def __hash__(self):
        return hash(self.to_tuple())

    def toggle(self, j):
        # Toggle bit at 1-based index j
        new_bf = self.bitfield ^ (1 << (j - 1))
        # Convert back to tuple for constructor
        return BitfieldRow([i for i in range(1, 64) if (new_bf & (1 << (i - 1))) != 0])

    def shifted(self, shiftup=0):
        # Shift all bits up by shiftup
        bf = self.bitfield << shiftup
        return BitfieldRow([i for i in range(1, 64) if (bf & (1 << (i - 1))) != 0])

    def __eq__(self, other):
        if isinstance(other, BitfieldRow):
            return self.to_tuple() == other.to_tuple()
        if isinstance(other, tuple):
            return self.to_tuple() == other
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, BitfieldRow):
            return self.to_tuple() < other.to_tuple()
        if isinstance(other, tuple):
            return self.to_tuple() < other
        return NotImplemented

    def __le__(self, other):
        if isinstance(other, BitfieldRow):
            return self.to_tuple() <= other.to_tuple()
        if isinstance(other, tuple):
            return self.to_tuple() <= other
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, BitfieldRow):
            return self.to_tuple() > other.to_tuple()
        if isinstance(other, tuple):
            return self.to_tuple() > other
        return NotImplemented

    def __ge__(self, other):
        if isinstance(other, BitfieldRow):
            return self.to_tuple() >= other.to_tuple()
        if isinstance(other, tuple):
            return self.to_tuple() >= other
        return NotImplemented

    def __ne__(self, other):
        if isinstance(other, BitfieldRow):
            return self.to_tuple() != other.to_tuple()
        if isinstance(other, tuple):
            return self.to_tuple() != other
        return NotImplemented
