class Tree(object):
    def __init__(self, name, children=None):
        self.name = name
        self.children = []
        if children is not None:
            assert(all(isinstance(c, Tree) for c in children))
            self.children = list(children)

    def __mul__(self, other):
        return Composite(self, other)

    def __pow__(self, p):
        return Power(self, p)


class Lagrange(Tree):
    def __init__(self, order, dimRange=1):
        Tree.__init__(self, "Lagrange")
        self.order = order
        self.dimRange = dimRange

    def __repr__(self):
        if self.dimRange == 1:
            return "Lagrange<" + str(self.order) + ">"
        else:
            return "Lagrange<" + str(self.order) + ">^" + str(self.dimRange)


class DG(Tree):
    def __init__(self, order, dimRange=1):
        Tree.__init__(self, "DG")
        self.order = order
        self.dimRange = dimRange

    def __repr__(self):
        if self.dimRange == 1:
            return "DG<" + str(self.order) + ">"
        else:
            return "DG<" + str(self.order) + ">^" + str(self.dimRange)


class Nedelec(Tree):
    def __init__(self, kind, order):
        Tree.__init__(self, "Nedelec")
        self.kind = kind
        self.order = order

    def __repr__(self):
        return "Nedelec<" + str(self.kind) + "," + str(self.order) + ">"


class RaviartThomas(Tree):
    def __init__(self, order):
        Tree.__init__(self, "RaviartThomas")
        self.order = order

    def __repr__(self):
        return "RaviartThomas<" + str(self.order) + ">"


class Composite(Tree):
    def __init__(self, *args, **kwargs):
        assert len(args) > 0
        Tree.__init__(self, "Composite", args)
        self.blocked = kwargs.get("blocked", True)
        self.layout = kwargs.get("layout", "lexicographic")

    def __repr__(self):
        return "(" + " * ".join(repr(c) for c in self.children) + ")"


class Power(Tree):
    def __init__(self, children, exponent, **kwargs):
        assert children is not None
        Tree.__init__(self, "Power", [children])
        assert len(self.children) == 1
        self.exponent = exponent
        self.blocked = kwargs.get("blocked", True)
        self.layout = kwargs.get("layout", "interleaved")

    def __repr__(self):
        if self.exponent == 1:
            return repr(self.children[0])
        else:
            return "[" + repr(self.children[0]) + "]^" + str(self.exponent)
