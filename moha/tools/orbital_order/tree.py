class Tree():
    def __init__(self, name='root', children=None):
        self.name = name
        self.parent = None
        self.children = []
        if children is not None:
            for child in children:
                self.add_child(child)
    def __repr__(self):
        return str(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def add_child(self, node):
        assert isinstance(node, Tree)
        self.children.append(node)
        node.parent = self

    def is_root(self):
        if self.parent is None:
            return True
        return False
    def is_leaf(self):
        if not self.children:
            return True
        return False
    def root(self):
        if self.is_root() is True:
            return self
        return self.parent.root()
    def dist_to_root(self):
        if self.is_root() is True:
            return 0
        return 1+self.parent.dist_to_root()
    def path_to_root(self):
        path = [self]
        if self.is_root() is True:
            return path
        path.extend(self.parent.path_to_root())
        return path
    def is_ancestor(self, node): # node is ancestor of self
        if node is self:
            return False
        elif node in self.path_to_root():
            return True
        return False
    def path_to_node(self, node):
        path1 = []
        temppath = self.path_to_root()
        for v in temppath:
            path1.append(v)
            if v is node:
                return (v, path1)
            if (node.is_ancestor(v) is True):
                mid = v
                break;
        temppath = node.path_to_root()
        path2 = []
        for v in temppath:
            path2.append(v)
            if v is mid:
                break
        path2.reverse()
        return (mid, path1+path2[1:])
