import unittest
import Attacks
import snap
import networkx as nx

class Test_Testing(unittest.TestCase):
    WBAtest = Attacks.WalkBasedAttack(0,1, "mygraph.txt")
    CBAtest = Attacks.CutBasedAttack("mygraph.txt")
    PAtest = Attacks.PassiveAttack("mygraph.txt")

    #tests for the walk-based attack

    #test convert subgraph function of SNAP
    def test_ConvertSubgraph(self):
        SetOfAccounts = [0,1,2,3]

        #create new graph G1
        G1 = snap.TUNGraph.New()

        #add 10 nodes to G1
        for i in range(10):
            G1.AddNode(i)
        
        #add a path from 0 to 9
        for i in range(9):
            G1.AddEdge(i, i+1)

        #create subgraph induced by nodes in SetOfAccounts
        TestH = G1.ConvertSubGraph(snap.TUNGraph, SetOfAccounts)
        nodes = []
        for NI in TestH.Nodes():
            nodes.append(NI.GetId())

        nodes.sort()
        self.assertEqual(SetOfAccounts,nodes)
        self.assertTrue(TestH.IsEdge(0,1))
        self.assertTrue(TestH.IsEdge(1,2))
        self.assertTrue(TestH.IsEdge(2,3))

    #test is edge function of SNAP
    def test_IsEdge(self):
        self.assertTrue(self.WBAtest.Graph.IsEdge(0,3))
        self.assertTrue(self.WBAtest.Graph.IsEdge(0,6))
        self.assertTrue(self.WBAtest.Graph.IsEdge(1,41))

     #test nodes() of snap
    def test_Nodes(self):
        Nodes = list(range(0,100))
        GNodes = []
        for NI in self.WBAtest.Graph.Nodes():
            GNodes.append(NI.GetId())

        GNodes.sort()
        self.assertEqual(Nodes, GNodes)

    #test the power set function of WBA
    def test_WBAPowerSet(self):
        values = [1,2,3]
        result = {(1,),(2,),(3,),(1,2),(1,3),(2,3),(1,2,3)}
        fvalues = set(self.WBAtest.powerset(values))
        self.assertEqual(result, fvalues)

    #test the power set function of PassiveAttack
    def test_PassivePowerSet(self):
        values = 3
        result = {(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)}
        fvalues = set(self.PAtest.powerset(values))
        self.assertEqual(result, fvalues)

    #test declareTargets function of WBA
    def test_WBA_declareTargets(self):
        amount = 4
        fvalues = self.WBAtest.declareTargets(4)
        Range = list(range(0,100))
        for i in range(len(fvalues)):
            self.assertIn(fvalues[i], Range)

    #test createUserAccounts function of WBA
    def test_WBA_createUserAccounts(self):
        amount = 4
        fvalues = self.WBAtest.createUserAccs(4)
        result = [100, 101, 102, 103]
        self.assertEqual(result, fvalues)

    #test checkExternalDegree function of WBA
    def test_WBA_checkExternalDegree(self):
        self.WBAtest.ExternalDegree = []
        self.WBAtest.ExternalDegree.append(10)
        self.WBAtest.ExternalDegree.append(10)
        self.WBAtest.ExternalDegree.append(10)
        self.WBAtest.ExternalDegree.append(0)

        self.assertTrue(self.WBAtest.checkExternalDegree([0,1,2]))
        self.assertFalse(self.WBAtest.checkExternalDegree([1,2,3]))

    #test constructEdges function of WBA
    def test_constructEdges(self):
        SetOfNodes = {100,101,102,103}
        target = 0 
        self.WBAtest.constructEdges(SetOfNodes, target)

        self.assertTrue(self.WBAtest.Graph.IsEdge(100,0))
        self.assertTrue(self.WBAtest.Graph.IsEdge(101,0))
        self.assertTrue(self.WBAtest.Graph.IsEdge(102,0))
        self.assertTrue(self.WBAtest.Graph.IsEdge(103,0))

    #test createInternalEdges function of WBA
    def test_createInternalEdges(self):
        SetOfAccounts = [100,101,102,103]
        self.WBAtest.TotalDegree = {}
        for i in SetOfAccounts:
            self.WBAtest.TotalDegree[i] = 0

        self.WBAtest.createInternalEdges(SetOfAccounts)
        #test if path from x_1 to x_4 exists
        self.assertTrue(self.WBAtest.Graph.IsEdge(100,101))
        self.assertTrue(self.WBAtest.Graph.IsEdge(101,102))
        self.assertTrue(self.WBAtest.Graph.IsEdge(102,103))

    #test getNeighborsWithDegree function of WBA
    def test_getNeighborsWithDegree(self):
        
        #in mygraph.txt node 3 has has 21 and 92 as a neighbor, which have degree 23
        node = 3
        degree = 23

        values = [21,92]
        fvalues = self.WBAtest.getNeighborsWithDegree(node, degree)
        
        self.assertEqual(values, fvalues)


    #tests for the cut based attack

    #test declareTargets function of WBA
    def test_CBA_declareTargets(self):
        amount = 4
        fvalues = self.CBAtest.declareTargets(4)
        Range = list(range(0,100))
        for i in range(len(fvalues)):
            self.assertIn(fvalues[i], Range)

    #test createUserAccs function of CBA
    def test_CBA_createUserAccs(self):
        values = [100,101,102,103]
        fvalues = self.CBAtest.createUserAccs(4)
        self.assertEqual(fvalues, values)

    #test convert subgraph function of Networkx
    def test_nx_subgraph(self):
        
        SetOfAccounts = [0,1,2,3]
        SetOfEdges = [(0,1), (1,2), (2,3)]
        #create new graph G1
        G = nx.Graph()
        
        #add nodes from 0 to 9, add a path from 0 to 9
        for i in range(9):
            G.add_edge(i, i+1)

        #create subgraph induced by nodes in SetOfAccounts
        H = nx.subgraph(G, SetOfAccounts)
        nodes = []
        edges = []
        for node in nx.nodes(H):
            nodes.append(node)

        for edge in nx.edges(H):
            edges.append(edge)

        nodes.sort()
        self.assertEqual(SetOfAccounts,nodes)
        self.assertEqual(SetOfEdges, edges)


    #test the gomory hu algorithm of networkx
    def test_gomory_hu(self):
        G = nx.Graph()
        G.add_edge(1,2)
        G[1][2]['capacity'] = 10

        G.add_edge(1,6)
        G[1][6]['capacity'] = 8

        G.add_edge(2,6)
        G[2][6]['capacity'] = 3

        G.add_edge(2,5)
        G[2][5]['capacity'] = 2

        G.add_edge(2,3)
        G[2][3]['capacity'] = 4

        G.add_edge(3,4)
        G[3][4]['capacity'] = 5

        G.add_edge(3,5)
        G[3][5]['capacity'] = 4

        G.add_edge(3,6)
        G[3][6]['capacity'] = 2

        G.add_edge(4,5)
        G[4][5]['capacity'] = 7

        G.add_edge(4,6)
        G[4][6]['capacity'] = 2

        G.add_edge(5,6)
        G[5][6]['capacity'] = 3

        #compute gomory hu tree of G
        GomoryHuTree = nx.gomory_hu_tree(G)

        resultG = nx.Graph()
        resultG.add_edge(1,2)
        resultG[1][2]['weight'] = 18

        resultG.add_edge(2,6)
        resultG[2][6]['weight'] = 17

        resultG.add_edge(6,5)
        resultG[6][5]['weight'] = 13

        resultG.add_edge(5,4)
        resultG[5][4]['weight'] = 14

        resultG.add_edge(5,3)
        resultG[5][4]['weight'] = 15

        TNodes = list(nx.nodes(GomoryHuTree)).sort()
        GNodes = list(nx.nodes(resultG)).sort()

        self.assertEqual(TNodes, GNodes)
        TEdges = set(nx.edges(GomoryHuTree))
        GEdges = set(nx.edges(resultG))
        self.assertEqual(TEdges,GEdges)

    #test the connected compoment function of networkx
    def test_connectec_components(self):
        G = nx.Graph()
        G.add_edge(1,2)
        G.add_edge(2,3)

        G.add_edge(4,5)
        G.add_edge(5,6)
        G.add_edge(6,4)

        values = [{1,2,3}, {4,5,6}]
        j = 0
        for i in nx.connected_components(G):
            self.assertEqual(values[j], i)
            j = j+1


    #test networkx is_isomorphic function
    def test_is_isomorphic(self):
        G = nx.Graph()
        H = nx.Graph()

        G.add_edge(1,2)
        G.add_edge(2,3)
        G.add_edge(3,4)
        G.add_edge(4,1)

        H.add_edge(5,6)
        H.add_edge(6,7)
        H.add_edge(7,8)
        H.add_edge(8,5)

        self.assertTrue(nx.is_isomorphic(G,H))

    
    #tests for the passive attack

    def test_compute_GofS_alpha(self):
        #potential coalition members
        CoalitionMembers = [300,301,302]
        self.PAtest.Graph.AddNode(300)
        self.PAtest.Graph.AddNode(301)
        self.PAtest.Graph.AddNode(302)

        #nodes connected to the potential coalition members
        self.PAtest.Graph.AddNode(401)
        self.PAtest.Graph.AddNode(402)
        self.PAtest.Graph.AddNode(403)
        self.PAtest.Graph.AddNode(404)
        self.PAtest.Graph.AddNode(405)
        self.PAtest.Graph.AddNode(406)
        self.PAtest.Graph.AddNode(407)
        self.PAtest.Graph.AddNode(408)
        self.PAtest.Graph.AddNode(409)

        #connecting the two sets of nodes
        self.PAtest.Graph.AddEdge(300,401)
        self.PAtest.Graph.AddEdge(300,402)
        self.PAtest.Graph.AddEdge(300,403)
        self.PAtest.Graph.AddEdge(300,404)

        self.PAtest.Graph.AddEdge(301,404)
        self.PAtest.Graph.AddEdge(301,405)
        self.PAtest.Graph.AddEdge(301,406)

        self.PAtest.Graph.AddEdge(302,405)
        self.PAtest.Graph.AddEdge(302,406)
        self.PAtest.Graph.AddEdge(302,407)
        self.PAtest.Graph.AddEdge(302,408)
        self.PAtest.Graph.AddEdge(302,409)

        GofS_alpha = self.PAtest.ComputeGofS_alpha(CoalitionMembers)
        values = {
            (0,) : 4,
            (1,) : 3,
            (2,) : 5,
            (0,1) : 6,
            (0,2) : 9,
            (1,2) : 6,
            (0,1,2) : 9
            }

        self.assertEqual(GofS_alpha, values)

if __name__ == '__main__':
    unittest.main()
