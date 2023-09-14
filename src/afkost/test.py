class _Kmer():
    def __init__(self):
        self.sequence = ""
        self.count = [0, 1]

    def __add__(self, item):
        if not isinstance(item, _Kmer):
            raise TypeError("`kmer` mist be a `_Kmer` instance")
        else:
            self.count[0] = self.count[0] + item.count[0]
            self.count[1] = self.count[1] + item.count[1]
            return self

    a = 1
    def return_count(self, value):
        print(self.count, value)
