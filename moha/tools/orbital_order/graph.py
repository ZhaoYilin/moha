class UGraph():
    '''
    Graph object for undirected graph

    :param graphlist: list of tuples, [(0,1),(2,3)], each pair stand for an edge with two vertex
    '''
    def __init__(self, graphlist):

        if not graphlist:
            self.is_weighted = False
        elif len(graphlist[0]) == 2:
            self.is_weighted = False
        elif len(graphlist[0]) == 3:
            self.is_weighted = True
        else:
            raise Exception('not a representation of a graph')

        self.vertex = set([graphlist[i][0] for i in range(len(graphlist))]
                          + [graphlist[i][1] for i in range(len(graphlist))])

        self.edge = []
        self.rep = []
        for e in graphlist:
            if e[0] < e[1]:
                self.edge.append((e[0], e[1]))
                if self.is_weighted is True:
                    self.rep.append((e[0],e[1],e[2]))
            else:
                self.edge.append((e[1], e[0]))
                if self.is_weighted is True:
                    self.rep.append((e[1],e[0],e[2]))

        self.adj = {}
        for v in self.vertex:
            self.adj[v] = []
            for edge in graphlist:
                if edge[0] == v:
                    if self.is_weighted is False:
                        self.adj[v].append(edge[1])
                    else:
                        self.adj[v].append((edge[1], edge[2]))
                if edge[1] == v:
                    if self.is_weighted is False:
                        self.adj[v].append(edge[0])
                    else:
                        self.adj[v].append((edge[0], edge[2]))



class Blossom():
    def __init__(self, rep): # rep: (base, loop path)
        self.base = rep[0]
        self.loop = rep[1]
        self.len = len(self.loop)
    def path_to_base(self, node):
        nid = self.loop.index(node)
        bid = self.loop.index(self.base)
        dist = nid - bid
        path = []
        if dist == 0:
            return [self.base]
        if dist %2 ==0 and dist<0:
            return self.loop[nid:bid+1]
        elif dist %2 ==0 and dist>0:
            for i in range(dist+1):
                path.append(self.loop[nid-i])
            return path
        elif dist %2 ==1 and dist<0:
            for i in range(nid+1):
                path.append(self.loop[nid-i])
            for i in range(self.len-bid):
                path.append(self.loop[-i-1])
            return path
        elif dist %2 ==1 and dist>0:
            for i in range(self.len-nid):
                path.append(self.loop[nid+i])
            for i in range(bid+1):
                path.append(self.loop[i])
            return path

