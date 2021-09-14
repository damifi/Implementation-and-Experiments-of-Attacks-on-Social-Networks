import math
import random
import array
import snap
import itertools
import time
import networkx as nx
import anytree as at


#Class for the walk-based attack
class WalkBasedAttack():
    
    #number of targets
    k = 0
    #list of every node's external degree
    ExternalDegree = None
    #list of every node's total degree
    TotalDegree = None
    #small integer constant, paper claims c=3 will suffice
    c = 0
    #number of targets
    b = 0
    #graph of the social network, used just for testing
    Graph = None
    #number of nodes of the social network
    SizeOfGraph = None
    #Integers, specifiying the bounds of the external degrees for every node
    DZero = None
    DOne = None
    #subgraph H, containing the user accounts the attacker created
    H = None
    #list of the created User accounts
    UserAccounts = None
    #Subsets that target users
    SubsetsTarget = None

    #constructor
    #input: 
    #DZero: lower bound for the external degree of every user account
    #DOne: upper bound for the external degree of every user account
    #path: path to the graph to be attacked
    def __init__(self, DZero: int, DOne: int, path: str) -> None:
        #load the edge list of the graph
        self.Graph = snap.LoadEdgeList(snap.TUNGraph, path, 0, 1)
        
        #Size of the graph
        self.SizeOfGraph = self.Graph.GetNodes()
        self.DOne = DOne
        self.DZero = DZero
    
    #Powerset computation
    #input: a list of IDs of nodes
    #output: a set containing the power set of all list elements (without the empty set) 
    #https://docs.python.org/3/library/itertools.html#itertools-recipes
    def powerset(self, ListOfNodes: list[int]):
        s = list(ListOfNodes)
        return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(1, len(s)+1))


    #create user accounts and insert them into the graph, currently at the end of the graph
    #returns a list of the ids of the new users
    #TODO use snaps function to give a node
    def createUserAccs(self, amount: int) -> list[int]:
        NewUserAccs = []
        #TODO: get the biggest ID and create accounts
        for i in range(self.SizeOfGraph, self.SizeOfGraph+amount):
            self.Graph.AddNode(i)
            NewUserAccs.append(i)

        self.UserAccounts = NewUserAccs
        return NewUserAccs

    #create a list of targets
    #for the tests I choose targets at random
    #input:
    #amount: Defines how many targets to choose
    #output:
    #list containing the IDs of the target
    def declareTargets(self, amount: int, ListOfTargets = None, chose = False) -> list[int]:
        if chose:
            #TODO: this can be removed
            return ListOfTargets
        else:
            return random.sample(range(0,self.SizeOfGraph), amount)

    #check for a set of nodes if every Node has extneral degree left. If it does the function returns True otherwise it returns False
    def checkExternalDegree(self, SetOfNodes: set[int]) -> bool:
        #iterate through all nodes
        for j in SetOfNodes:
            #if current node has no external degree left return false
            if(self.ExternalDegree[j] <= 0):
                return False
        return True

    # creates links between a set of nodes and a target
    #input:
    #Nodes: set of node IDs
    #target: ID of the target
    def constructEdges(self, Nodes: set[int], target: int) -> None:
        #iterate through all nodes in the list and construct edges to the target
        for i in Nodes:
            self.Graph.AddEdge(i,target)

    # Function to create edges from H to G-H
    def fillEdges(self, SetOfAccounts: list[int], SetOfTargets: list[int]) -> None:
        #list that contains keys of the dictionary to delete, this is the case if an user account is uniquely connected
        #to a target, but it also has external degree left. This means the user account must be connected to other nodes in G
        #meaning the target cannot be uniquely identified. 
        KeysToDelete = []
        for Key in self.SubsetsTarget.keys():
            if len(Key) == 1:
                if self.ExternalDegree[Key[0]] > 0:
                    KeysToDelete.append(Key)
        #Delete every node from the dictionary that has a target that cannot be uniquely identified        
        for Key in KeysToDelete:
            del self.SubsetsTarget[Key]

        #connect node i to random accounts until it has Delta_i many neighbors in G
        for i in SetOfAccounts:
            Sample = random.sample(range(0, self.SizeOfGraph-len(SetOfAccounts)), self.ExternalDegree[i])
            # if sample contains a target compute another sample, chance should be very small in large social networks
            while not set(Sample).isdisjoint(SetOfTargets):
                Sample = random.sample(range(0, self.SizeOfGraph-len(SetOfAccounts)), self.ExternalDegree[i])
            for j in Sample:
                self.Graph.AddEdge(i, j)
            self.ExternalDegree[i] = 0

    # Function to create the interal edges in H to increase the chance of recovery
    def createInternalEdges(self, SetOfAccounts: list[int]) -> None:
        #create the path from x1 to xk
        i = 0
        #list containing nodes from the path x1 to xk
        PathList = []
        while(i < len(SetOfAccounts)-1):
            #add an edge between x_i and x_i+1
            self.Graph.AddEdge(SetOfAccounts[i], SetOfAccounts[i+1])
            #increase the total degree of x_i and x_i+1
            self.TotalDegree[SetOfAccounts[i]] = self.TotalDegree[SetOfAccounts[i]] + 1
            self.TotalDegree[SetOfAccounts[i+1]] = self.TotalDegree[SetOfAccounts[i+1]] + 1
            PathList.append((SetOfAccounts[i], SetOfAccounts[i+1]))
            i = i + 1

        #create other edges(xi,xj)
        #create all combinations of users
        for i in itertools.combinations(SetOfAccounts,2):
            if i not in PathList:
                if(random.randint(0,1) == 0):
                    self.Graph.AddEdge(i[0], i[1])
                    self.TotalDegree[i[0]] = self.TotalDegree[i[0]] + 1
                    self.TotalDegree[i[1]] = self.TotalDegree[i[1]] + 1

    #Creates the subgraph H
    #input:
    #SetOfAccounts: list of node IDs, corresponding to the user accounts the attacker created
    #SetOfTargets: list of node IDs, corresponding to the targets the attacker has chosen
    #c: maximum size of the sets N_j
    def createH(self, SetOfAccounts: list[int], SetOfTargets: list[int], c = 3) -> int: 
        StartA = time.process_time()
        StartB = time.time()
        print("Starting to construct subgraph H!")
        self.c = c
        #List of sets of all possible subsets
        PowerSetOfAccounts = list(self.powerset(SetOfAccounts))

        #remove all subsets that are bigger than c
        SetOfDistinctSets = [s for s in PowerSetOfAccounts if len(s) <= self.c]

        #the algorithm takes the reverse order of the Setofdistinctsets to choose the subsets that attack the target
        SetOfDistinctSets.reverse()

        #Dict that contains the subsets that connect to the targets, N_j attacks target j 
        self.SubsetsTarget = {}


        # 1) Compute the external degree for every node:
        self.ExternalDegree = {}
        for i in SetOfAccounts:
            self.ExternalDegree[i] = int(random.uniform(self.DZero, self.DOne))

        #Create a list containing the total degree of every node, start at external degree
        self.TotalDegree = self.ExternalDegree.copy()
        
        # 2) Choose the subsets that are connected to the targets 
        i = 0
        while i < len(SetOfTargets):
            # stop the algorithm if no more subsets are available
            if not SetOfDistinctSets:
               break

            # check if the current subset has external degree left
            if(self.checkExternalDegree(SetOfDistinctSets[0])):
               for j in SetOfDistinctSets[0]:
                   self.ExternalDegree[j] = self.ExternalDegree[j] - 1

               #construct the edges between the set of user accounts and the target
               self.constructEdges(SetOfDistinctSets[0], SetOfTargets[i])
               self.SubsetsTarget[(SetOfDistinctSets[0])] = SetOfTargets[i]
               #Remove the used set
               SetOfDistinctSets.pop(0)
               i += 1
            else:
               SetOfDistinctSets.pop(0)
                
        
        
        # 3) Create further links from H to G
        self.fillEdges(SetOfAccounts, SetOfTargets)

        # 4) Create the internal edges in H such that H is easier to recover
        self.createInternalEdges(SetOfAccounts)

        #Get subgraph H
        self.H = self.Graph.ConvertSubGraph(snap.TUNGraph, SetOfAccounts)

        print(f"Construction of Subgraph H completed. CPU running time: {time.process_time() - StartA}. Wall-clock running time: {time.time() - StartB}")
        return len(self.SubsetsTarget)


    # Function to find neighbors of a node with the given degree. Returns a list containing the neighbors
    #input:
    #Node: ID of the node to check for neighbors
    #Degree: degree of neighbors we are looking for
    #output:
    #list of all neighbors of Node with the correct degree
    def getNeighborsWithDegree(self, Node: int, Degree: int) -> list[int]:
        # v contains nodes that have the correct amount of neighbours, i.e have the correct TotalDegree
        v = []
        # iterate through the neighbors of the given Node and put all neighbors with the correct degree into list v
        for Neighbor in self.Graph.GetNI(Node).GetOutEdges():
            node = self.Graph.GetNI(Neighbor)
            if node.GetOutDeg() == Degree:
                v.append(Neighbor)
        
        return v

    #Function to perform the structure test of subgraph H
    # input: 
    # Node: nodeID as integer that is the potential next node on the path
    # Path: a path to the current node
    # output: Boolean indicating, whether the edges back on the path exist 
    def checkEdges(self, Node: int, Path: list[int]) -> bool:

        #remove the root from the path
        Path = Path[1:]
        l = len(Path)
        
        #check for duplicates, a path containing duplicates cannot correspond to the path we are looking for,
        #as all user accounts are unique
        PathDub = Path
        PathDub.append(Node)
        if len(PathDub) != len(set(PathDub)):
            return False
        
        #iterate through all nodes in the path and perform the structure test:
        #edge (alpha_i, alpha_l) exists if and only if edge (x_i, x_l) exists
        #return true if all edges are correct 
        for i in range(l):
            if(self.H.IsEdge(self.UserAccounts[i], self.UserAccounts[l])) == (self.Graph.IsEdge(Path[i], Node)):
                pass
            else:
                return False
        return True


    #Function that perfroms the recovery algorithm of the walk-based attack
    def recoverH(self):

        StartA = time.process_time()
        StartB = time.time()

        print("\nStarting to recover H!\n")
        #attack cannot succeed for k <= 1, therefore return fail
        if( len(self.UserAccounts) <= 1):
            return (False, time.process_time() - StartA, 0, 0)

        #Iterator to go through TotalDegrees
        Iterator = iter(self.TotalDegree)
        Degree = self.TotalDegree[next(Iterator)]
        NextDegree = self.TotalDegree[next(Iterator)]

        #use tree from anytree as the rooted search tree
        from anytree import Node, RenderTree, PreOrderIter
        
        #Create dummy node alpha*
        Root = Node("Dummy")

        #First step: create the first leafs in T
        #iterate over all nodes and find nodes with degree x_1
        for node in self.Graph.Nodes():
            #check if current node has the desired degree
            if node.GetOutDeg() == Degree:

                #adding nodes with TotalDegree of x_1 to the tree
                TempNode = Node(node.GetId(), parent = Root)                

        #Tree size after all potential first nodes are added, used in experiments to make a guess how large the tree typically is
        TreeSizeFirstIteration = len(at.findall(Root))

        #make a list of total Degrees
        DegreeTest = []
        for i in self.TotalDegree.values():
                DegreeTest.append(i)

        #counter of how many unique k-paths exist
        UniquePaths = 0
        #Leaf node of the first length-k path found
        LeafNode = None
        #iterating through all leaf nodes
        for node in PreOrderIter(Root):
            if node.is_leaf:
                #if the leafs path has length k the path has reached its maximum length, therefore continue with
                #the next potential paths
                if node.depth == len(DegreeTest):
                    #check how many unique potential paths have been found
                    #if more than one are found, the attack did not succeed
                    UniquePaths = UniquePaths + 1
                    if UniquePaths >= 2:
                        TreeSize = len(at.findall(Root))
                        return (False, time.process_time() - StartA, TreeSize, TreeSizeFirstIteration)
                    #keep track of the first leaf node that 
                    LeafNode = node
                    continue
                
                #NextDegree is the next degree we are looking for, indicated by tree depth
                NextDegree = DegreeTest[node.depth]

                # Get all neighbors of current leaf node with the correct degree
                V = self.getNeighborsWithDegree(node.name, NextDegree)

                # List that contains the nodes of the path to the current leaf node
                PathToLeaf = []
                for n in node.path:
                    PathToLeaf.append(n.name)

                #check if the edges (f(alpha_i), v) is an edge iff (x_i, x_l+1) is an edge
                for NB in V:
                    if self.checkEdges(NB, PathToLeaf):
                        temp = Node(NB, parent = node)
            
        #H is a list of nodes that correspond to the nodes in subgraph H
        H = []
        if UniquePaths == 1:
            print("Attack succeeded. Found a unique path in G: ")
            for n in LeafNode.path:
                H.append(n.name)

            if H[1:] == self.UserAccounts:

                print(f"Recovery of Subgraph H completed. CPU running time: {time.process_time() - StartA}. Wall-clock running time: {time.time() - StartB}\n\n")
                TreeSize = len(at.findall(Root))
                return (True, time.process_time() - StartA, TreeSize, TreeSizeFirstIteration)
        else:
            print("Was not able to find a unique path. Attack did not work.")
            print(f"Recovery of Subgraph H completed. CPU running time: {time.process_time() - StartA}. Wall-clock running time: {time.time() - StartB}\n\n")
            TreeSize = len(at.findall(Root))
            return (False, time.process_time() - StartA, TreeSize, TreeSizeFirstIteration)
        
        print(f"Something went wrong. Wall-clock running time: {time.time() - StartB}\n\n")


class CutBasedAttack():
    #Number of users we wish to target
    b = None
    #List of Users we want to attack
    W = None
    #number of new user accounts
    k = None
    #List of new user accounts
    X = None
    #minimum degree in H
    DeltaOfH = None
    #minimum cut in H
    GammaOfH = None
    #Graph of the social network
    Graph = None
    #Size of the social network
    SizeOfGraph = None
    #Subgraph H
    SubgraphH = None
    #Networkx Graph of the social network
    NetworkXGraph = None

    #constructor
    #input:
    #path: path to the graph to perform the attack on
    def __init__(self, path: str) -> None:
        self.NetworkXGraph = nx.read_edgelist(path, nodetype = int)
        self.SizeOfGraph = self.NetworkXGraph.order()


    #create a list of targets, currently the nodes are chosen at random
    #input:
    #amount: int to specify how many targets to choose
    def declareTargets(self, amount: int) -> list[int]:
        return random.sample(range(0,self.SizeOfGraph), amount)

    #create user accounts and insert them into the graph
    #input:
    #amount: int to specify how many accounts to create
    #returns a list of the ids of the new users
    def createUserAccs(self, amount: int) -> list[int]:
        NewUserAccs = []
        for i in range(self.SizeOfGraph, self.SizeOfGraph+amount):
            self.NetworkXGraph.add_node(i)
            NewUserAccs.append(i)
        return NewUserAccs

    #function to create internal edges of H with probability 1/2
    def createInternalEdges(self) -> None:
        #all possible pairs of the user accounts:
        for i in list(itertools.combinations(self.X, 2)):
            #if random.choice([0,1]) == 1:
            if random.random() <= 0.50:
                self.NetworkXGraph.add_edge(i[0], i[1])

    #chooses users in X that will be attached to targets
    def chooseUsers(self) -> None:
        #choose b nodes in X arbitrarily
        users = random.sample(self.X, self.b)
        for i in range(len(self.W)):
            print(f"Adding edge from attacker {users[i]} to target {self.W[i]} ")
            self.NetworkXGraph.add_edge(users[i], self.W[i])


    #Function to create the subgraph H and insert it into G
    #input:
    #NumberOfTargets: int to specify how many targets to attack
    def constructH(self, NumberOfTargets: int) -> None:

        StartA = time.process_time()
        StartB = time.time()

        # Step 1)
        #initialize all parameters
        self.b = NumberOfTargets
        self.W = self.declareTargets(self.b)
        self.k = 3* self.b + 3
        self.X = self.createUserAccs(self.k)
        print(f"b: {self.b} \n W: {self.W} \n k: {self.k} \n X: {self.X}")
        self.createInternalEdges()

        # Step 2)
        # choose users to attach to the targets
        self.chooseUsers()
        #Create subgraph H
        self.SubgraphH = nx.subgraph(self.NetworkXGraph, self.X)


        print(f"Construction of Subgraph H completed. CPU running time: {time.process_time() - StartA}. Wall-clock running time: {time.time() - StartB}\n\n")


    #Function that implements the recovery algorith of the cut-based attack
    def recoverH(self):                                

        StartA = time.process_time()
        StartB = time.time()

        # 1)
        # create a graph in the format of networkx
        #Construct Gomory-Hu tree of G
        #set the edge weight of each edge to 1
        nx.set_edge_attributes(self.NetworkXGraph, 1, "capacity")
        print("Computing Gomory-Hu tree...")
        from networkx.algorithms import flow
        #compute the Gomory-Hu tree of G

        T = nx.gomory_hu_tree(self.NetworkXGraph, flow_func = flow.dinitz, capacity='capacity')
        print("Gomory-Hu tree computation complete.")

        # 2)
        #delete all edges with weight at most b from T, producing forest T'
        EdgesToRemove = [(a,b) for a,b, attributes in T.edges(data=True) if attributes['weight'] <= self.b]
        T.remove_edges_from(EdgesToRemove)
                
        #iterate through all components of size k -> Sets S
        S = []
        for i in nx.connected_components(T):
            #print(i)
            if len(i) == self.k:
                S.append(self.NetworkXGraph.subgraph(i).copy())

        #if no S can be found, the attack cannot succeed.
        if not len(S):
            print("No Subset in forest T found. Attack was not successful.")
            print(f"Recovery of Subgraph H failed. CPU running time: {time.process_time() - StartA}. Wall-clock running time: {time.time() - StartB}\n\n")
            return (0, time.process_time() - StartA)

        #test whether Sets S are isomorpic to H
        print("Searching for isomorphic Graphs...")
        IsomorphCounter = 0
        FoundNodes = None
        for s in S:
            if nx.is_isomorphic(self.SubgraphH, s):
                IsomorphCounter = IsomorphCounter + 1
                FoundNodes = s.nodes()

        print(f"Found {IsomorphCounter} isomorphic subsets to H.")

        #if only one isomorph component has been found the attack succeeded.
        if(IsomorphCounter == 1):
            print("Attack succeeded. Found a single subgraph in G that is isomorph to H")
            print(f"List of created accounts: {self.X}")
            print(f"Nodes : {FoundNodes}")
            print(f"Recovery of Subgraph H succeeded. CPU running time: {time.process_time() - StartA}. Wall-clock running time: {time.time() - StartB}\n\n")
            return (1, time.process_time() - StartA)
        else:
            print("Attack failed.")
            print(f"Recovery of Subgraph H failed. CPU running time: {time.process_time() - StartA}. Wall-clock running time: {time.time() - StartB}\n\n")
            return (0, time.process_time() - StartA)

class PassiveAttack:
    
    #Graph of the social network
    Graph = None
    #Set of the coalition members
    X = None
    #number of coalition members
    k = None
    #Subgraph induced by nodes of the coaltion members
    SubgraphH = None
    #Coalition member that recruits other coalition members
    StartingNode = None
    #Total degree of each coalition member
    TotalDegree = None
    #Dictionary that contains the number of who the coalition members are connected to
    GofS = {}

    #Constructor 
    #input:
    #path: path to the graph's edge list
    def __init__(self, path: str) -> None:
        #load the graph form the specified path
        self.Graph = snap.LoadEdgeList(snap.TUNGraph, path, 0, 1)
        self.X = []

    #function to select the starting node, who recruits all other nodes
    #input:
    #ID: ID of the first node
    def selectNode(self, ID: int) -> None:
        self.X.append(ID)

    #recruit all other nodes, return 0 if amount neighbors are recruited.
    def selectNeighbors(self, Amount: int) -> int:
        counter = 0
        
        #iterate through neighbours of the first coalition member and add them to the coalition
        for Neighbor in self.Graph.GetNI(self.StartingNode).GetOutEdges():
            if counter >= Amount:
                return 0
            else:
                self.X.append(Neighbor)
                counter = counter + 1


    #Return all neighbours of the first node x_1 and the correct degree we are looking for
    #input:
    #DegreeIndex: index of the degree we are looking for
    #FirstNode: ID of the first node
    #output: List of all neighbours with the correct degree
    def FindNeighbors(self, DegreeIndex: int, FirstNode: int) -> list[int]:
        v = []
        #always look for neighbors of x_1 as x_1 is the only node connected to all members
        #and get all neighbours with correct next degree indicated by degree index
        for Neighbor in self.Graph.GetNI(FirstNode).GetOutEdges():
            node = self.Graph.GetNI(Neighbor)
            if node.GetOutDeg() == self.TotalDegree[DegreeIndex]:
                v.append(Neighbor)

        return v



    #powerset
    #https://docs.python.org/3/library/itertools.html#itertools-recipes
    def powerset(self, Length: int) -> set[int]:
        s = range(Length)
        return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(1, len(s)+1))

    #checks whether the nodeID given as paramater and current sequence of nodes has the same structure as the one in H 
    def checkEdges(self, Node: int, Sequence: list[int]) -> bool:
       
        #remove the dummy node for the sequence
        Sequence = Sequence[1:]
        l = len(Sequence)

        #for each node in the sequence and the potential node to be added to the sequence,
        #perform the structure test
        for i in range(l):
            if(self.SubgraphH.IsEdge(self.X[i], self.X[l])) == (self.Graph.IsEdge(Sequence[i], Node)):
                pass
            else:
                return False
        return True
        


    #Compute the amount of unique neighbors that are connected to the nodes in Nodes and stores these numbers in a dictionary
    #input: list of nodes to compute g(S) for
    #output: g(S) for the coalition members
    def computeGofS(self, Nodes: list[int]) -> None:
        #powerset of all possible indices
        S = list(self.powerset(len(Nodes)))
        #Dict containing the neighbor sets for all given nodes
        NeighborDict = {}
        #go through all nodes
        for Node in Nodes:
            #create a new set for the neighbors
            Neighbors = set()
            #go through all neighbors of the current node
            #for NB in self.NetworkxGraph.neighbors(Node):
            for NB in self.Graph.GetNI(Node).GetOutEdges():
                #add the neighbor to the set if it is outside of H
                if NB not in Nodes:
                    Neighbors.add(NB)
            NeighborDict[Node] = Neighbors

        #go through all sets in S
        for Tuple in S:
            
            TempSet = set()
            #go through all IDs in Tuple
            for ID in Tuple:
                #Get the union of all sets indicated by the ID
                TempSet = TempSet.union(NeighborDict[Nodes[ID]])
            #The number of users the set of coalition members is connected to
            self.GofS[Tuple] = len(TempSet)

    #similar to ComputeGofS, this function computes the amount of neighbors of nodes in Nodes and returns these numbers in a dictionary
    #input: list of nodes to compute g(S)_alpha for
    #output: g(S)_alpha for the current sequence in T
    def ComputeGofS_alpha(self, Nodes: list[int]):
        # all possible subsets from 0...k-1
        S = list(self.powerset(len(Nodes)))

        NeighborDict = {}
        for Node in Nodes:
            #create a new set for the neighbors
            Neighbors = set()
            #go through all neighbors of the current node
            for NB in self.Graph.GetNI(Node).GetOutEdges():
                #add the neighbor to the set if it is outside of H
                if NB not in Nodes:
                    Neighbors.add(NB)
            NeighborDict[Node] = Neighbors

        #go through all sets in S
        GofS_alpha = {}
        for Tuple in S:
            
            TempSet = set()
            #go through all IDs in Tuple
            for ID in Tuple:
                #Get the union of all sets indicated by the ID
                TempSet = TempSet.union(NeighborDict[Nodes[ID]])
            #The number of users the set of coalition members is connected to
            GofS_alpha[Tuple] = len(TempSet)
        return GofS_alpha

    #Function to perform the passive attack
    #input:
    #CoalitionSize: number of the nodes that form the coalition
    def attack(self, CoalitionSize : int):
        print("Starting the passive attack")
        StartA = time.process_time()
        StartB = time.time()

        self.k = CoalitionSize
        #Get candidates for possible x1 (StartingNode), i.e nodes with degree >= NumberToRecruit
        Candidates = []
       
        #go through all nodes and create a list of all nodes that have enough neighbours to be x_1
        for NI in self.Graph.Nodes():
            if NI.GetOutDeg() >= (self.k - 1):
                Candidates.append(NI.GetId())

        #Select a starting node, currently at random
        self.StartingNode = random.choice(Candidates)
        self.selectNode(self.StartingNode)

        #Select all the other neighbors that form a coalition
        self.selectNeighbors(CoalitionSize-1)

        #Compute the total degree of the coalition members 
        self.TotalDegree = []

        #Get the total degree of every coalition member
        for Node in self.X:
            self.TotalDegree.append(self.Graph.GetNI(Node).GetOutDeg())
        

        #Save the subgraph H as an new graph
        self.SubgraphH = self.Graph.ConvertSubGraph(snap.TUNGraph, self.X)
      
        print(f"Construction of Subgraph H completed. CPU running time: {time.process_time() - StartA}. Wall-clock running time: {time.time() - StartB}\n\n")
        
        #computing the optimization from the paper g(S)
        self.computeGofS(self.X)

        print("Starting Recovery of Subgraph H")

        #using anytree as a tree for the rooted search tree
        from anytree import Node, RenderTree, PreOrderIter, findall
        
        #initialize T and create dummy node alpha*
        Root = Node("Dummy")

        #looking for nodes with the same degree as x_1 and creating the first nodes in the rooted search tree
        for node in self.Graph.Nodes():
            if node.GetOutDeg() == self.TotalDegree[0]:
                TempNode = Node(node.GetId(), parent = Root)
            

        #variable containing the results of the attack
        AttackResults = []
        #First value of results is the size of the tree after the first iteration
        AttackResults.append(len(at.findall(Root)))


        #variable to keep track of the amount of unique sequences
        #UniqueSequences = 0
        #go through all leaf nodes
        for node in PreOrderIter(Root):
            if node.is_leaf:
                #if current leaf is at depth k continue with the next leaf
                if node.depth == self.k:
                    continue

                #For the current leaf node get all the neighbors with the correct degree
                Neighbors = self.FindNeighbors(node.depth, node.path[1].name)

                # List that contains the path to the current leaf node
                SequenceToLeaf = []

                for n in node.path:
                    SequenceToLeaf.append(n.name)

                #check if the edges (f(alpha_i), v) is an edge iff (x_i, x_l+1) is an edge
                #and if so add the neighbour as a leaf to the tree
                for NB in Neighbors:
                    if self.checkEdges(NB, SequenceToLeaf):
                        temp = Node(NB, parent = node)

        
        #Second value of the tree is the size of the tree after all sequences are constructed
        AttackResults.append(len(at.findall(Root)))

        #Counter that counts how many sequences were found that match subgraph H
        CorrectSequences = 0
        #List of the nodes that correspond to the nodes of the coalition, if the attack was successful
        CorrectSequence = []
        #Go through all leaf nodes
        for Node in Root.leaves:
            #Sequence that is to be checked
            Sequence = []
            #Check if the current node's depth is equal to the coalition size
            if Node.depth == self.k:
                #go through all nodes in the current leaf's path
                for node in Node.path:
                    #add the nodes to the sequence list
                    Sequence.append(node.name)
                #remove root from the sequence
                Sequence.pop(0)
                #Compute GofS_alpha, this is the refined algorithm from the paper
                GofS_alpha = self.ComputeGofS_alpha(Sequence)
                #If GofS from the original nodes is equal to GofS_alpha of the current Sequence this must mean
                #with high probability the current sequence is a match to be the subgraph H
                if self.GofS == GofS_alpha:
                    #Count how many sequences with correct GofS_alpha are found
                    CorrectSequences = CorrectSequences + 1
                    #Store the potential match of H in CorrectSequence
                    CorrectSequence = Sequence
                else:
                    #print("GofS != GofS_alpha")
                    pass
        #print(f"Number of sequences that match H: {CorrectSequences}")
        if CorrectSequences == 1:
            print(f"CorrectSequence: {CorrectSequence}, Coalition: {self.X}")
            #print(list(self.GofS.values()))
            #print(list(self.GofS.values()).count(1))
            # 3rd value is if the attack succeeded
            AttackResults.append(True)

            # 4th value is the amount of nodes which can be uniquely identified, that is, if a node is connected to a subset
            # of colition members, and g(S) = 1
            AttackResults.append(list(self.GofS.values()).count(1))
            
            
        else:
            print("Attack did not succeed.")
            # 3rd value: Attack failed
            AttackResults.append(False)
            # 4th value, number ob targets = 0
            AttackResults.append(0)


        print(f"Recovery of Subgraph H completed. CPU running time: {time.process_time() - StartA}. Wall-clock running time: {time.time() - StartB}\n\n")
        # last value is the time the algorith took
        AttackResults.append(time.process_time() - StartA)
        return AttackResults



        