from collections import defaultdict


class Require:
    _called = defaultdict(lambda: False)

    @property
    def called(self):
        return self._called

    def __init__(self, requirements):
        self.requirements = requirements

    def __call__(self, f):
        for func in self.requirements:
            if not self.called[func]:
                f.

        def wrap(*args, **kwargs):
            return f(*args, **kwargs)
        self.called[f] = True
        return wrap


class Example:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def play(self):
        self.a, self.b = self.b, self.a

    @Require([play])
    def test(self):
        print(self.a, self.b)


if __name__ == '__main__':
    example = Example(5, 10)
    example.test()
