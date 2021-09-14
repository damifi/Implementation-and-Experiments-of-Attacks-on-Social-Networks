import snap

class Utilities:

    Graph = None


    #path to a graph in .txt format
    def __init__(self, path: str):
        print("reading graph")
        self.Graph = snap.LoadEdgeList(snap.TNGraph, path, 0, 1)
        print("Done loading the graph")
       
       


    #function to convert the graph into an undirected one
    #nodes are renamed to 0,1,...,n-1
    #input:
    #name: name the graph should be saved as
    def ConvertToUndirected(self, name: str) -> None:
        UndirectedGraph = self.Graph.ConvertGraph(snap.TUNGraph, True)
        UndirectedGraph.SaveEdgeList(name)


    #function to print the number of nodes and the degree sequence of G,
    #prints the average degree and the amount of edges per node
    def GetAverageDegree(self):

        totaldegree = 0
        for item in self.Graph.GetDegCnt():
            print("%d nodes with degree %d" % (item.GetVal2(), item.GetVal1()))
            totaldegree = totaldegree + (item.GetVal2() * item.GetVal1())

        print(f"totaldegree: {totaldegree}")
        print(f"averagedegree: {totaldegree/self.Graph.GetNodes()}")
        print(f"edges per node: {self.Graph.GetEdges()/self.Graph.GetNodes()}")