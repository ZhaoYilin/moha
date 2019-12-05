from moha.tools.orbital_order.graph import UGraph
from copy import deepcopy
import logging
logger = logging.getLogger('blossom')


def edge_reg(e):  # make sure the left endpoint is smaller
    if e[0] < e[1]:
        return e
    return (e[1], e[0])


def list_shift(l, n):  # shift the list in a loop fashion to right, the original input list is changed!
    for i in range(n):
        l.insert(0, l.pop(-1))
    return l


def reduce_graph(graph, btuple):
    new_edges = []
    blossom_ele = list(btuple)
    for e in graph.edge:
        ee = list(e)
        if (e[0] in blossom_ele) and (e[1] not in blossom_ele):
            ee[0] = btuple[0]
            if edge_reg(tuple(ee)) not in new_edges:
                new_edges.append(edge_reg(tuple(ee)))
        elif (e[1] in blossom_ele) and (e[0] not in blossom_ele):
            ee[1] = btuple[0]
            if edge_reg(tuple(ee)) not in new_edges:
                new_edges.append(edge_reg(tuple(ee)))
        elif (e[1] in blossom_ele) and (e[0] in blossom_ele):
            pass
        else:
            if e not in new_edges:
                new_edges.append(e)
    return UGraph(new_edges)


class DressedGraph:
    def __init__(self, graph):
        self.graph = graph
        self.n = len(graph.vertex)  # number of vertex of original graph

        if self.n == 0:
            logger.warning('Empty graph')
            maxweight = 0
        else:
            maxweight = max(graph.rep, key=lambda x: x[2])[2]
        self.m = len(graph.edge)  # number of edge of original graph

        self.edge = graph.edge
        self.weight = {(e[0], e[1]): e[2] for e in graph.rep}  # edge weight correspondence

        self.z = {i: maxweight / 2 for i in graph.vertex}  # valid for all vertex v and blossom in any level
        # blossom is indexed by the tuple rep, which is stored in self.status
        self.status = {i: [(i,)] for i in graph.vertex}
        self.contain = {i: [i] for i in
                        graph.vertex}  # valid for top v, for v not in top level, its elements are for diff levels
        # by status together with contain, we can know the over all structure of the nested blossoms
        # status: {1: [(1,),(1,2),(3,1)], 2: [(2,),(1,2),(3,1)], 3:[(3,),(3,1)]}, each item of dict track the blossoms from bottom to top
        # contain: {1:[1,2,3,4], 2:[3,4]...} each item of dict track all basic vertex contained in some high level effective vertex

        self.match = UGraph([])  # relative matching on top level
        self.rgraph = UGraph(self.graph.edge)  # reduced graph
        # quantities below clear every phase
        self.stf = {i: 'F' for i in graph.vertex}  # valid for all real v, no matter level it is
        self.stfb = {i: None for i in graph.vertex}  # valid for top blossom, labeled as v, None stands no blossom in v
        self.parent = {i: None for i in graph.vertex}  # valid for top level, b labeled as v
        self.mark = {i: False for i in graph.edge}  # defined on edge of top level, relative mark
        self.checkq = []  # check queue

        self.static_status = None # desinged for keep blossom info for the final all-expansion

    def is_edge(self, edge):
        e = edge_reg(edge)
        if e in self.edge:
            return True
        return False

    def blossoms_cover_edge(self, edge): # given edge based on original vertex, return list of blossom tuples the edge in
        e = edge_reg(edge)
        blist = []
        i = -1
        while(1):
            if self.static_status[e[0]][i] == self.static_status[e[1]][i]:
                blist.append(self.static_status[e[0]][i])
                i -= 1
            else:
                break
        return blist


    def is_equal(self, edge): # check top level edge: whether they are equal for z and w
        e = edge_reg(edge)

        # ? how to choose the edge weight for the merging edge
        ee_cand = self.edge_toold(e)

        for ee in ee_cand:
            z = 0
            if self.static_status is not None:
                z = 0
                blist = self.blossoms_cover_edge(ee)
                for b in blist:
                    z += self.z[b]
                #     logger.debug('blossom including %s, z is %s:' %(b,self.z[b]))#
                # logger.debug('z for two end is: %s, %s, edge(%s,%s) weight is %s' % (self.z[ee[0]], self.z[ee[1]],
                #                                                                     ee[0],ee[1],self.weight[ee]))  #
            if self.weight[ee] == z + self.z[ee[0]] + self.z[ee[1]]:
                return True
        return False

    def top_vertice(self):  # set of labels in top level as vertex
        topset = set()
        for v in self.status.items():
            topset.add(v[1][-1][0])
        return topset

    def blossom_level(self, v):  # 0 for bare nodes on top level, n for n nested blossom for basic v vertex
        return len(self.status[v]) - 1

    def top_rep(self, v): # given a basic vertex, return the top level rep of the blossom containing v
        return self.status[v][-1][0]

    def edge_tonew(self, edge):  # from initial into graph of top level
        e = edge_reg(edge)
        res = edge_reg((self.top_rep(e[0]), self.top_rep(e[1])))
        if res[0] == res[1]:
            return None  # original edges within blossom
        return res

    def edge_toold(self, edge): # from edge in top level to original edge, a list of edge
        e = edge_reg(edge)
        result = []
        for u in self.contain[e[0]]:
            for v in self.contain[e[1]]:
                newe = edge_reg((u, v))
                if (self.is_edge(newe) is True) and (newe not in result):
                    result.append(newe)
        return result

    def stf_check(self, v):  # v in top level: is it S or T
        if self.blossom_level(v) == 0:
            return self.stf[v]
        else:
            return self.stfb[v]

    def stf_set(self, v, label):  # v in top level
        logger.debug('set top level vertex %s with label %s' % (v, label))
        if self.blossom_level(v) == 0:
            self.stf[v] = label
            self.stfb[v] = None
        else:
            self.stfb[v] = label
            for child in self.contain[v]:
                self.stf[child] = label

    def find_root(self, v):
        if self.parent[v] is None:
            return v
        return self.find_root(self.parent[v])

    def path_to_root(self, v):
        path = [v]
        if self.parent[v] is None:
            return path
        path.extend(self.path_to_root(self.parent[v]))
        return path

    def path_to_node(self, v, w): # if v w in one tree, return path to form blossom, with base in the first element
        path1 = []
        pathv = self.path_to_root(v)
        pathw = self.path_to_root(w)
        mid = pathv[-1]
        for node in pathv:
            path1.append(node)
            if node in pathw:
                mid = node
                break
        path2 = []
        for v in pathw:
            path2.append(v)
            if v is mid:
                break
        path2.reverse()
        path = path1 + path2[1:]
        midid = path.index(mid)
        path = list_shift(path, len(path) - midid)
        return tuple(path)

    def compose_blossom(self, btuple):  # the first one is the rep element of the blossom, which is the based when added
        logger.debug('compose new blossom made use of top level: %s' % list(btuple))
        assert set(btuple).issubset(self.top_vertice())

        # adjust mark edges
        logger.debug('rgraph edge: %s' % self.rgraph.edge)
        logger.debug('the mark dict: %s'%self.mark)
        for e in self.rgraph.edge:
            if (e[0] in btuple) and (e[1] not in btuple):
                # boolv = self.mark[edge_reg(e)]
                del self.mark[edge_reg(e)]
                self.mark[edge_reg((e[1], btuple[0]))] = False
            elif (e[1] in btuple) and (e[0] not in btuple):
                # boolv = self.mark[edge_reg(e)]
                del self.mark[edge_reg(e)]
                self.mark[edge_reg((e[0], btuple[0]))] = False
            elif (e[1] in btuple) and (e[0] in btuple):
                del self.mark[edge_reg(e)]
            else:
                pass  # edges outsisde the new blossom has no effect




        self.rgraph = reduce_graph(self.rgraph, btuple)
        logger.debug('new top level graph generated, with vertices %s'%self.rgraph.vertex)
        self.match = reduce_graph(self.match, btuple)

        for bov in btuple[1:]:  # every element in btuple is itself a vertex or blossom
            self.contain[btuple[0]].extend(self.contain[bov])
        for v in self.contain[btuple[0]]:
            self.status[v].append(btuple)
        self.z[btuple] = 0  # add a new blossom into z dict which keep all the vertex and blossom

        for v in btuple:
            if self.blossom_level(v) != 0:  # v is already a blossom before merging
                self.stfb[v] = None
        self.stf_set(btuple[0], 'S')  # always in S position when a new blossom is generated
        for v in btuple:  # clear items in check queue who are now not in top level
            try:
                self.checkq.remove(v)
            except ValueError:
                pass
        self.checkq.append(btuple[0])

        # update parent outside the blossom
        for v in self.top_vertice():
            if self.parent[v] in btuple[1:]:
                self.parent[v] = btuple[0]

        logger.debug('parent after composing blossom is: %s'%self.parent)
        logger.debug('ending of composing new blossom: %s' % list(btuple))
        logger.debug('contain of top level: %s'%self.contain)
        logger.debug('status: %s' % self.status)

    def expand_blossom(self, vrep, delz = True):  # v is the rep of the blossom at top level
                                                # delz is False for final expand-all
        if self.blossom_level(vrep) == 0:  # compatible with vertex expansion
            return None
        btuple = self.status[vrep][-1]
        logger.debug('expand blossom with vertex contained: %s'%list(btuple))
        # check whether there is a child attached to this blossom
        child = None
        for v in self.top_vertice():
            if vrep == self.parent[v]:
                child = v
                break

        logger.debug('z: %s'%self.z)#
        # logger.debug('the status:%s'% self.status)
        if delz is True:
            del self.z[btuple]
        for v in self.contain[btuple[0]]:
            self.status[v].pop(-1)
        for v in btuple[1:]:
            for u in self.contain[v]:
                self.contain[vrep].remove(u)

        logger.debug('the status:%s'%self.status)
        logger.debug('the contain:%s' % self.contain)
        self.stfb[vrep] = None

        self.graph_update()  # expand the graph

        # update parent and stf(b)
        father = self.parent[vrep]
        for u in btuple:
            self.parent[u] = None
            self.stf_set(u, 'F')
        if father is not None:
            for ve in self.rgraph.adj[father]:
                if ve in btuple and self.is_equal(edge_reg((ve, father))): # or (delz is False)
                    start = ve
                    break
            if child is None:

                self.parent[start] = father
                self.stf_set(start, 'T')
            else:

                for ve in self.rgraph.adj[child]:
                    if ve in btuple and self.is_equal(edge_reg((ve, child))) : #or (delz is False)
                        end = ve
                        break
                self.parent[child] = end
                # sid = btuple.index(start)
                eid = btuple.index(end)
                blist = list(btuple)

                blist = list_shift(blist, len(blist) - eid)  # eid at index 0
                # logger.debug('the  blist: %s' % blist)  #
                sid = blist.index(start)
                eid = blist.index(end)
                if (eid - sid) % 2 == 1:
                    blist.reverse()
                    blist = list_shift(blist, 1)
                    # logger.debug('the blist: %s' % blist)  #
                logger.debug('the start, end, and the blist: %s,%s,%s' % (start, end, blist))  #
                for j, ve in enumerate(blist):

                    if ve == start:
                        self.parent[ve] = father
                        logger.debug('set %s parent as %s'%(ve, father))
                        self.stf_set(ve, 'T')
                        break
                    else:
                        self.parent[ve] = blist[j + 1]
                        logger.debug('set %s parent as %s' % (ve, blist[j+1]))
                        if j % 2 == 0:
                            self.stf_set(ve, 'T')
                        else:
                            self.stf_set(ve, 'S')

        logger.debug('outermost stf of blossoms: %s'%self.stfb)
        logger.debug('parent state after blossom exposion: %s'%self.parent)

        # update mark
        new_mark = {}
        for e in self.rgraph.edge:
            if e in self.mark and (e[0] not in btuple) and (e[1] not in btuple):
                new_mark[e] = self.mark[e]
            else:
                if (self.parent[e[0]] == e[1]) or (self.parent[e[1]] == e[0]):
                    new_mark[e] = True
                else:
                    new_mark[e] = False
        self.mark = new_mark
        # logger.debug('mark state after blossom expansion: %s'%self.mark)

        new_edges = self.match.edge
        # update relative match
        # check is there is a match line:
        blist = list(btuple)
        if vrep in self.match.vertex:
            outer = self.match.adj[vrep][0]
            for member in btuple:
                #logger.debug('member %s, outer %s, in edges %s, equal %s'%(member, outer, edge_reg((member, outer)) in self.rgraph.edge
                                                                         #  ,self.is_equal(edge_reg((member, outer)))))#
                if (edge_reg((member, outer)) in self.rgraph.edge):
                    # if delz is False:
                    #     logger.debug('member %s, outer %s, equal %s, blossom included %s'
                    #                  %(member, outer,self.is_equal(edge_reg((member, outer))), self.blossoms_cover_edge(edge_reg((member, outer))) ))#
                    if (self.is_equal(edge_reg((member, outer)))): #or (not delz): # delz is in last phase, there is nonzero z
                        # blossom to be expanded and hence is_euqal doesn't work
                        inner = member
                        if self.parent[outer] == member:
                            break
            # logger.debug('remove edge (%s,%s) from new match'%(vrep, outer))
            new_edges.remove( edge_reg( (vrep, outer) ) )
            # logger.debug('add edge (%s,%s) into new match' % (inner, outer))
            new_edges.append(edge_reg((inner, outer)))
            # choose match edges within blossom based on matched (inner,outer) edge

            rem = btuple.index(inner)
            # logger.debug('inner point and its index is: %s, %s' %(inner, rem))
            blist = list_shift(blist, len(btuple) - rem)
            # logger.debug('blossom loop from base point: %s'%blist)
        for i, _ in enumerate(blist):
            if i % 2 == 1:
                new_edges.append(edge_reg((blist[i], blist[i + 1])))
        self.match = UGraph(new_edges)
        logger.debug('ending of expansion blossom: %s' % list(btuple))

    def graph_update(self):  # reconstruct the rgraph based on current self.status
        graph_edges = set()
        for e in self.graph.edge:
            newe = self.edge_tonew(e)
            if newe is not None:
                graph_edges.add(newe)
        self.rgraph = UGraph(list(graph_edges))

    def phase_start(self): # initialize and reset some states for a new phase
        self.stf = {i: 'F' for i in self.graph.vertex}  # valid for all real v, no matter level it is
        self.stfb = {i: None for i in
                     self.graph.vertex}  # valid for top blossom, labeled as v, None stands no blossom in than v
        for v in self.rgraph.vertex:
            if self.blossom_level(v) > 0:
                self.stfb[v] = 'F'
        self.parent = {i: None for i in self.graph.vertex}
        # valid for top level, b labeled as v, None stands no parent at now Or not in the top level
        self.mark = {e: False for e in self.rgraph.edge}
        self.checkq = []

    def return_aug_path(self, v, w):
        logger.debug('find aug path connect two trees connection is %s:%s' % (v, w))
        p1 = self.path_to_root(v)
        p2 = self.path_to_root(w)
        p1.reverse()
        return p1 + p2  # aug path from v to w as a list

    def find_aug_path(self):
        for e in self.match.edge:
            if e not in self.rgraph.edge:
                raise Exception('illegal match')
        exposed_vertices = self.rgraph.vertex - self.match.vertex
        for v in exposed_vertices:
            if self.stf_check(v) != 'S':
                self.stf_set(v, 'S')
                self.checkq.append(v)

        for e in self.match.edge:
            self.mark[e] = True


        while (len(self.checkq) != 0):
            logger.debug('checkq now: %s' % self.checkq)
            v = self.checkq.pop(0)
            logger.debug('try develop tree from vertex: %s' % v)
            # logger.debug('the adj matrix is %s' % self.rgraph.adj)
            if v not in self.rgraph.vertex:# deal with graph with separated points
                pass
            else:
                for w in self.rgraph.adj[v]:

                    e = edge_reg((v, w))
                    # logger.debug('edge %s-%s, mark %s euqal %s'%(v,w,self.mark[e],self.is_equal(e)))
                    if self.mark[e] is False and (self.is_equal(e) is True):
                        logger.debug('try edge (%s,%s)' % (v, w))
                        if self.stf_check(w) == 'S':  # w is already in the forest with S position
                            if self.find_root(v) != self.find_root(w):
                                path = self.return_aug_path(v, w)  # tbd
                                return path
                            else:
                                btuple = self.path_to_node(v, w)
                                self.compose_blossom(btuple)
                                break
                        elif self.stf_check(w) == 'F':
                            self.add_to_forest(v, w)
                        self.mark[e] = True

        return []

    def add_to_forest(self, v, w):  # insert w and its partner pairs to v
        x = self.match.adj[w][0]
        self.parent[w] = v
        self.parent[x] = w
        self.stf_set(w, 'T')
        self.stf_set(x, 'S')
        self.checkq.append(x)

    def phase(self):
        logger.info('matches now are: %s' % self.match.edge)
        logger.info('-----------------------')
        logger.info('a new phase start...')

        self.phase_start()
        while (1):
            p = self.find_aug_path()
            if p:
                self.alt_path(p)
                return False
            logger.info('a new subphase for weight adjustment start...')
            t = self.reweight()
            logger.info('finished adjusment type: %s' % t)
            if t == 1 or t == 0:
                return True

    def phases(self):
        if self.n == 0:
            return False
        while (1):
            r = self.phase()

            # expand all S blossoms with zero z after each stage?

            zerosb = []
            for v in self.rgraph.vertex:
                if (self.stfb[v] == 'S') and self.z[self.status[v][-1]]==0:
                    zerosb.append(v)
            while (zerosb):
                b = zerosb.pop(0)
                btuple = self.status[b][-1]
                self.expand_blossom(b)
                for v in btuple:
                    if (self.stfb[v] == 'S') and self.z[self.status[v][-1]]==0:
                        zerosb.append(v)

            logger.info('end of one phase round')


            if r is True:
                # expand all blossoms at very last of all?
                logger.info('expand all blossoms finally')
                self.static_status = deepcopy(self.status)
                logger.debug('final static status: %s'%self.static_status)
                bqueue = []
                for v in self.top_vertice():
                    if self.stfb[v] is not None:
                        bqueue.append(v)
                while(bqueue):
                    logger.debug('bqueue:%s'%bqueue)#
                    b = bqueue.pop(0)
                    btuple = self.status[b][-1]
                    self.expand_blossom(b, False)
                    for v in btuple:
                        if self.stfb[v] is not None:
                            bqueue.append(v)
                # self.check_optimum()
                return True

    def reweight(self): # reweight all z, if no possible path can be found
        # 1
        d1 = None
        v1 = []
        for v in self.graph.vertex:
            if self.stf[v] == 'S':
                if d1 is None:
                    d1 = self.z[v]
                    v1 = [v]
                elif d1 > self.z[v]:
                    d1 = self.z[v]
                    v1 = [v]
                elif d1 == self.z[v]:
                    v1.append(v)
        # 2
        d2 = None
        e2 = []
        for e in self.graph.edge:
            # logger.debug('d2 review: edge %s(%s)-%s(%s), value is %s' % (e[0],self.stf[e[0]], e[1],self.stf[e[1]],  self.weight[e]))  #
            if {self.stf[e[0]], self.stf[e[1]]} == {'F', 'S'}:
                tem = self.z[e[0]] + self.z[e[1]] - self.weight[e]

                if d2 is None:
                    d2 = tem
                    e2 = [e]
                elif tem < d2:
                    d2 = tem
                    e2 = [e]
                elif tem == d2:
                    e2.append(e)

        # 3
        d3 = None
        e3 = []
        for e in self.graph.edge:
            if (self.stf[e[0]] == 'S') and (self.stf[e[1]] == 'S') and (self.top_rep(e[0]) != self.top_rep(e[1])):
                tem = (self.z[e[0]] + self.z[e[1]] - self.weight[e]) / 2
                if d3 is None:
                    d3 = tem
                    e3 = [e]
                elif tem < d3:
                    d3 = tem
                    e3 = [e]
                elif tem == d3:
                    e3.append(e)

        # 4
        d4 = None
        b4 = []
        for v in self.graph.vertex:
            if self.stfb[v] == 'T':
                tem = self.z[self.status[v][-1]] / 2
                if d4 is None:
                    d4 = tem
                    b4 = [v]
                elif tem < d4:
                    d4 = tem
                    b4 = [v]
                elif tem == d4:
                    b4.append(v)

        l = [(x, i + 1) for i, x in enumerate([d1, d2, d3, d4]) if x is not None]
        if not l:
            return 0  # all four types of delta are not defined: no free vertex
        logger.debug('the candidates of four delta: %s'%l)
        delta, d = min(l, key=lambda x: x[0])  # 1,2,3,4 which is the minimum
        logger.debug('the ajust value is %s' % delta)

        for v in self.graph.vertex:
            if self.stf[v] == 'S':
                self.z[v] = self.z[v] - delta
            elif self.stf[v] == 'T':
                self.z[v] = self.z[v] + delta
        for v in self.top_vertice():
            if self.stfb[v] == 'S':
                self.z[self.status[v][-1]] = self.z[self.status[v][-1]] + 2 * delta
            elif self.stfb[v] == 'T':
                self.z[self.status[v][-1]] = self.z[self.status[v][-1]] - 2 * delta

        logger.debug('the z value after adjustment: %s' % self.z)

        if d == 2:

            for ee in e2:

                e = self.edge_tonew(ee)
                logger.debug('new equal edge: (%s,%s)' %e)
                for i in (0, 1):
                    if self.stf_check(e[i]) == 'S':
                        self.checkq.append(e[i])
        elif d == 3:

            for ee in e3:
                e = self.edge_tonew(ee)
                logger.debug('new equal edge: (%s,%s)' %e)
                for i in (0, 1):
                    if self.stf_check(e[i]) == 'S':
                        self.checkq.append(e[i])
        elif d == 4:
            logger.debug('zero z T blossom: %s' % b4)
            logger.debug('parent states: %s'%self.parent)
            for b in b4:
                self.expand_blossom(b)
            for v in self.top_vertice():
                if self.stf_check(v) == 'S':
                    self.checkq.append(v)

        logger.debug('stf state: %s'%self.stf)


        return d

    def alt_path(self, path):
        logger.info('the augmented path to be alternated: %s' % path)
        new_edges = self.match.edge
        for i in range(len(path) - 1):
            if i % 2 == 0:
                new_edges.append((path[i], path[i + 1]))
            else:
                new_edges.remove(edge_reg((path[i], path[i + 1])))
        self.match = UGraph(new_edges)

    def z_on_edge(self, edge): # given original edge, return the z value of the edge, by considering all nonzero blossom
        z = 0
        if self.static_status is not None:
            blist = self.blossoms_cover_edge(edge)
            for b in blist:
                z += self.z[b]
        z += self.z[edge[0]]+self.z[edge[1]]
        return z


    def check_optimum(self): # final result check based on theorem
        for e in self.graph.edge:
            assert self.z_on_edge(e) >= self.weight[e]
        for e in self.match.edge:
            assert self.z_on_edge(e) == self.weight[e]
        free_vertices = self.graph.vertex - self.match.vertex
        for v in free_vertices:
            assert self.z[v] == 0
        # check the final result satisfing all conditions of z



def find_maximum_matching(graph_edges, debug = False):
    '''
    find match with maximum total weights on undirected weighted graph

    :param graph_edges: edge list in the form : [(2,3,6),(5,4,3),(4,2,7)], in each tuple, the first two elements are edge
        end vertex, and the last term is the weight of corresponding edge
    :param debug: boolean, default False, if set as True, the results would be check based on the theorem
    :returns: UGraph object, self.edge attribute give the mathcing edge tuple list
    '''
    st = DressedGraph(UGraph(graph_edges))
    st.phases()
    if debug is True:
        st.check_optimum()
    return st.match


def find_perfect_maximum_matching(graph_edges):
    '''
    find perfect matching with maximum total weights on undirected weighted graph

    :param graph_edges: edge list in the form : [(2,3,6),(5,4,3),(4,2,7)], in each tuple, the first two elements are edge
       end vertex, and the last term is the weight of corresponding edge
    :param debug: boolean, default False, if set as True, the results would be check based on the theorem
    :returns: UGraph object, self.edge attribute give the mathcing edge tuple list
    '''
    weights_sum = sum(e[2] for e in graph_edges)
    reweight_graph = []
    for e in graph_edges:
        reweight_graph.append((e[0],e[1],e[2]+weights_sum))
    st = DressedGraph(UGraph(reweight_graph))
    st.phases()
    return st.match
