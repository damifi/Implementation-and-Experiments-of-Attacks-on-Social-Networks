# Implementation-and-Experiments-of-Attacks-on-Social-Networks

Installation:
1. Install SNAP for python: python -m pip install snap-stanford
(https://snap.stanford.edu/snappy/index.html)
2. Install NetWorkx for python: pip install networkx
(https://pypi.org/project/networkx/)
3. Install Anytree for python: pip install anytree
(https://anytree.readthedocs.io/en/latest/installation.html)

Usage:
main.py is the entry point of the program from which the attacks are executed.
To execute from command line: python main.py

From now on the user interacts with the console:
First, the user can choose what attack to perform.
Then, a graph has to be chosen (for example: facebook_combined.txt).
Finally, some parameters can be chosen, depending on the attack.

The results are written to .txt files in the same directory, after the attacks are executed.

For testing:
mygraph.txt is required to be in the same folder as Testing.py to perform the tests
(python -m unittest -v Testing.py)
