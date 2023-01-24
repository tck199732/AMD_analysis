class InfiniteDict(dict):
    def __getitem__(self, key):
        return self.get(key) if key in self else self.setdefault(key, InfiniteDict())     