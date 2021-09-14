import Attacks
import snap
import time
import Utilities



#ask what attack to perform
choice = int(input("What attack to perform? 0 = Walk-based attack, 1 = Cut-based attack, 2 = Passive attack"))


#walk-based attack
if choice == 0:
    #Ask the user to specify a path to a graph.
    PathToGraph = input("Enter a path to a graph for the attacks: ")
    HowManyRuns = int(input("How many runs of the attack to perform?"))
    dzero = int(input("dzero?"))
    done = int(input("done?"))
    c = int(input("c?"))

    Start = time.process_time()
    Starta = time.time()

    for i in range(1,11):
        NumberOfAccounts = i
        #WalkBased Attack
        #counter of successful runs
        SuccessCounter = 0
        TotalTargets = 0
        TotalRecoveryTime = 0
        TreeSize = 0
        TreeSizeFirstIteration = 0
        successfulTreeSize = 0
        f = open("WBA_"+PathToGraph+"_Accounts_"+str(NumberOfAccounts)+"_Dzero_"+str(dzero)+"_Done_"+str(done)+".txt", "a")
        for i in range(HowManyRuns):
            Success = None
            WBA = Attacks.WalkBasedAttack(dzero, done, PathToGraph)
            SetOfTargets = WBA.declareTargets(300)
            SetOfAccounts = WBA.createUserAccs(NumberOfAccounts)
            HowManyTargetsCanBeAttacked = WBA.createH(SetOfAccounts, SetOfTargets, c)
            Results = WBA.recoverH()
            if Results[0]:
                Success = "Run " +str(i)+ " was succesful. Number of targets that can be attacked: " +str(HowManyTargetsCanBeAttacked)+", size of the tree: "+str(Results[2])+"\n"
                SuccessCounter = SuccessCounter + 1
                TotalTargets = TotalTargets + HowManyTargetsCanBeAttacked
                TotalRecoveryTime = TotalRecoveryTime + Results[1]
                TreeSize = TreeSize + Results[2]
                TreeSizeFirstIteration = TreeSizeFirstIteration + Results[3]
                successfulTreeSize = successfulTreeSize + Results[2]
            else:
                Success = "Run " +str(i)+ "failed." + " Size of the tree: "+str(Results[2])+"\n"
                TotalRecoveryTime = TotalRecoveryTime + Results[1]
                TreeSize = TreeSize + Results[2]
                TreeSizeFirstIteration = TreeSizeFirstIteration + Results[3]
            f.write(Success)

        SuccessRate = SuccessCounter / HowManyRuns * 100
        f.write("Number of accounts created: " + str(NumberOfAccounts))
        f.write("\nDZero: " + str(dzero))
        f.write("\nDOne: " + str(done))
        f.write("\nc: " + str(c))
        f.write("\nGraph that was attacked: " + PathToGraph)
        f.write("\nSuccessful attacks: " + str(SuccessCounter))
        f.write("\nSuccess rate: " + str(SuccessRate))
        if SuccessCounter > 0:
            f.write("\nAverage amount of attacked targets: "+str(TotalTargets/SuccessCounter))
            f.write("\nAverage size of the tree for only successful runs: "+str(successfulTreeSize/SuccessCounter))
        else:
            f.write("No successful runs. No targets could be attacked.")
            f.write("No successful runs. average size of tree for successful runs undifined.")

        f.write("\nAverage size of the tree: "+str(TreeSize/HowManyRuns))
    
        f.write("\nAverage size of the tree after the first iteration: " +str(TreeSizeFirstIteration/HowManyRuns))
        f.write("\nTotal recovery time: " +str(TotalRecoveryTime))
        f.write("\nAverage time the recovery algoritm took: " +str(TotalRecoveryTime/HowManyRuns))
        f.close()


 #CutBased attack
if choice == 1:
   
    #Ask the user to specify a path to a graph.
    PathToGraph = input("Enter a path to a graph for the attacks: ")
    HowManyRuns = int(input("How many runs of the attack to perform?"))
    NumberOfTargets = int(input("How many targets to attack?"))

    Start = time.process_time()
    Starta = time.time()
    CBSuccesses = 0
    TotalRecoveryTime = 0

    f = open("CBA_"+PathToGraph+"_Targets_"+str(NumberOfTargets)+".txt", "a")
    for i in range(HowManyRuns):
        CBA = Attacks.CutBasedAttack(PathToGraph)
        CBA.constructH(NumberOfTargets)
        result = CBA.recoverH()
        TotalRecoveryTime = TotalRecoveryTime + result[1]
        CBSuccesses = CBSuccesses + result[0]

    f.write("\nGraph that was attacked: " + PathToGraph)
    f.write("\nSuccessful attacks: " + str(CBSuccesses))
    f.write("\nNumber of performed attacks: " +str(HowManyRuns))
    f.write("\nSuccess rate: " + str(CBSuccesses/HowManyRuns * 100))
    f.write("\nTotalRecoveryTime " + str(TotalRecoveryTime))
    f.write("\nAverage recovery time: " +str(TotalRecoveryTime/HowManyRuns))
    f.close()

#Passive attack
if choice == 2:
    PathToGraph = input("Enter a path to a graph for the attacks: ")
    HowManyRuns = int(input("How many runs of the attack to perform?"))

    Start = time.process_time()
    Starta = time.time()

    for i in range(1,11):
        SuccessCounter = 0
        TotalTargets = 0
        TreeSizeFirstIteration = 0
        TreeSize = 0
        successfulTreeSize = 0
        TotalRecoveryTime = 0
        f = open("Passive_"+PathToGraph+"_Accounts_"+str(i)+".txt", "a")
        for run in range(HowManyRuns):
            Success = None
            #Passive attack
            PA = Attacks.PassiveAttack(PathToGraph)
            Results = PA.attack(i)

            if Results[2]:
                Success = "Run " +str(run)+ " was succesful. Number of targets attacked: " +str(Results[3])+", size of the tree: "+str(Results[1])+"\n"
                SuccessCounter = SuccessCounter + 1
                TreeSizeFirstIteration = TreeSizeFirstIteration + Results[0]
                TreeSize = TreeSize + Results[1]
                TotalTargets = TotalTargets + Results[3]
                TotalRecoveryTime = TotalRecoveryTime + Results[4]
                successfulTreeSize = successfulTreeSize + Results[1]

            else:
                Success = "Run " +str(run)+ "failed." + " Size of the tree: "+str(Results[1])+"\n"
                TotalRecoveryTime = TotalRecoveryTime + Results[4]
                TreeSize = TreeSize + Results[1]
                TreeSizeFirstIteration = TreeSizeFirstIteration + Results[0]
            f.write(Success)


        SuccessRate = SuccessCounter / HowManyRuns * 100
        f.write("\nGraph that was attacked: " + PathToGraph)
        f.write("\nSuccessful attacks: " + str(SuccessCounter))
        f.write("\nSuccess rate: " + str(SuccessRate))
        if SuccessCounter > 0:
            f.write("\nAverage amount of attacked targets: "+str(TotalTargets/SuccessCounter))
            f.write("\nAverage size of the tree for only successful runs: "+str(successfulTreeSize/SuccessCounter))
        else:
            f.write("No successful runs. No targets could be attacked.")
            f.write("No successful runs. average size of tree for successful runs undifined.")

        f.write("\nAverage size of the tree: "+str(TreeSize/HowManyRuns))
    
        f.write("\nAverage size of the tree after the first iteration: " +str(TreeSizeFirstIteration/HowManyRuns))
        f.write("\nTotal recovery time: " +str(TotalRecoveryTime))
        f.write("\nAverage time the recovery algoritm took: " +str(TotalRecoveryTime/HowManyRuns))
        f.close()

else:
    print("No valid input.")




print(f"Total CPU running time: {time.process_time()-Start} seconds")
print(f"Total wall-clock running time: {time.time() - Starta} seconds")



