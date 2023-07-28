#define _CRT_SECURE_NO_WARNINGS

#include<iostream>
#include<fstream>
#include<string>
#include<cstring>
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<set>
#include<unordered_map>
#include<utility>
#include<functional>
#include<cstdio>
#include<ctime>
#include<chrono>
#include<algorithm>
#include<cmath>
#include<queue>

using namespace std;

//GLOBALS
double epsilon; //Epsilon Threshold Value
double phi; //Phi Threshold Value
int maxFoundSize; //Size of Largest Found Structure
unordered_map<int, int> mpvMap; //Maps Maximum Positive Values
unsigned maxDeg; //Max Degree of nodes in G

/***
* Writes a string to file
*/
void writeToFile(string s, char* filename) {
	ofstream writeFile;
	writeFile.open(filename, ios::app);
	if (writeFile.is_open()) {
		writeFile << s << "\n";
		writeFile.flush();
	} 
	writeFile.close();
}

/***
* Writes a string to file
*/
void writeToFile(string s, string filename) {
	ofstream writeFile;
	writeFile.open(filename, ios::app);
	if (writeFile.is_open()) {
		writeFile << s << "\n";
		writeFile.flush();
	} 
	writeFile.close();
}

//CLASSES
/***
* The Node Class
*/
class Node {
	public:
		int id;
		vector<int> posNeighbours; //Positive Edges
		vector<int> negNeighbours; //Negative Edges
	
		vector<int> activePosNeighbours; //Positive Edges in the Search Space
		vector<int> activeNegNeighbours; //Negative Edges in the Search Space
		vector<Node*> activeNeighbours; //Active Neighbours in the Search Space
	
		int positiveCoreNumber; //PCN
		int corePosEdges; //Positive Number for PosCore Decomp
		int coreTotEdges; //Total Number for PosCore Decomp
        int negativeLimit; //Negative Limit value
        int negativePotential; //Negative Potential value
        bool active; //Is the node still active
	
		//Functions for inserting neighbour node information
		void insertNeighbour(int node, bool sign) {
			if (sign) {
				posNeighbours.push_back(node);
			} else {
				negNeighbours.push_back(node);	
			}
		}
	
		void insertPosNeighbour(int node) {
			posNeighbours.push_back(node);
		}
	
		void insertNegNeighbour(int node) {
			negNeighbours.push_back(node);	
		}
    
        //Calculate the Negative Potential
        void calculateNP() {
            double lowerVal = floor((double) activePosNeighbours.size() * (1.0 - epsilon));
            negativeLimit = (int) lowerVal;
            negativePotential = min((int) activeNegNeighbours.size(), negativeLimit);
        }
};

/***
* Helps to keep track of the epsilon values while performing BK
*/
class CandidateNode {
	public:
		int id;
		int posEdges;
        int negativeLimit;
	
        //Calculates Epsilon of current state
		double getEpsilon(int totEdges) {
			return (double) posEdges / (double) totEdges;
		}
    
        //Checks if the node exceeds the NL
        bool checkLimit(int size) {
            if ((size - posEdges) > negativeLimit) {
                return false;
            } else {
                return true;
            }
        }
};

/***
* Helps to keep track of the phi value while performing BK
*/
class CandidateClique {
	public:
		unordered_map<int, CandidateNode> nodes;
		int negEdges;
		int totEdges;
        int currentNP;
	
		//Inserts the Initial Query Node
		void insertNode(CandidateNode qNode) {
			nodes.insert(make_pair(qNode.id, qNode));
		}	
    
        void insertNode(Node node) {
            CandidateNode newCNode;
			newCNode.id = node.id;
            newCNode.posEdges = 0;
			nodes.insert(make_pair(newCNode.id, newCNode));
            currentNP += node.negativePotential;
        }
		
        double getPhi() {
			return (double) negEdges / (double) totEdges;
		}
    
        int size() {
			return nodes.size();
		}
    
        //Checks if the current P cup R violates NP
        bool checkTotalNP(int newNP) {
            double target = ceil((double) totEdges * phi);
            if ((newNP + currentNP) < (int) target) {
                return false;
            } else {
                return true;
            }
        } 
    
		//Inserts a new node and adjusts the relevent counters
		void adjustCounters(Node newNode) {
			int posEdgeChange = 0;
			int negEdgeChange = 0;
			unordered_map<int, CandidateNode>::iterator cIt;
            
            
            for (auto a : nodes) {
                cIt = nodes.find(a.second.id);
                if (binary_search(newNode.activePosNeighbours.begin(), newNode.activePosNeighbours.end(), a.first)) {
                    posEdgeChange++;
                    cIt->second.posEdges++;
                } else {
                    negEdgeChange++;
                }
            }
            negEdges = negEdges + negEdgeChange;
			totEdges = totEdges + negEdgeChange + posEdgeChange;
			CandidateNode newCNode;
			newCNode.id = newNode.id;
			newCNode.posEdges = posEdgeChange;
			nodes.insert(make_pair(newCNode.id, newCNode));
        }
    		
        //Inserts a new node, adjusts the counters and checks if the counters
        //are all correct s.t. the current state satisfies requirements
        bool adjustAndCheckCounters(Node newNode, int gamma, bool &limitBreaker) {
            
            if (size() + 1 < gamma) {
                adjustCounters(newNode);
                return false;
            }
            
            bool passable = true;
            int posEdgeChange = 0;
            int negEdgeChange = 0;
            int individualSize = size();
            unordered_map<int, CandidateNode>::iterator cIt;
            
            for (auto a : nodes) {
                cIt = nodes.find(a.second.id);
                if (binary_search(newNode.activePosNeighbours.begin(), newNode.activePosNeighbours.end(), a.first)) {
                    posEdgeChange++;
                    cIt->second.posEdges++;
                } else {
                    negEdgeChange++;
                }
                if (!cIt->second.checkLimit(individualSize)) {
                    limitBreaker = false;
                    return false;
                }
                if (passable && cIt->second.getEpsilon(individualSize) < epsilon) {
                    passable = false;
                }
            }
            
            negEdges = negEdges + negEdgeChange;
			totEdges = totEdges + negEdgeChange + posEdgeChange;
			CandidateNode newCNode;
			newCNode.id = newNode.id;
			newCNode.posEdges = posEdgeChange;
            newCNode.negativeLimit = newNode.negativeLimit;
            if (!newCNode.checkLimit(individualSize)) {
                limitBreaker = false;
                return false;
            }
            if (passable && newCNode.getEpsilon(individualSize) < epsilon) {
                passable = false;
            }
			nodes.insert(make_pair(newCNode.id, newCNode));
            currentNP += newNode.negativePotential;
            
            if (passable && getPhi() < phi) {
                passable = false;
            }
            
            return passable;
        }
    
};

//COMPARITORS
struct NodeComp {
    bool operator()(const Node &n1, const Node &n2) {
        if (n1.corePosEdges == n2.corePosEdges) {
            if (n1.coreTotEdges == n2.coreTotEdges) {
                return n1.id < n2.id;
            }
            return (n1.coreTotEdges < n2.coreTotEdges);
        }
        return (n1.corePosEdges < n2.corePosEdges);
    }
    
    bool operator()(const Node* &n1, const Node* &n2) {
        if ((*n1).corePosEdges == (*n2).corePosEdges) {
            if ((*n1).coreTotEdges == (*n2).coreTotEdges) {
                return (*n1).id < (*n2).id;
            }
            return ((*n1).coreTotEdges < (*n2).coreTotEdges);
        }
        return ((*n1).corePosEdges < (*n2).corePosEdges);
    }
    
    bool operator()(const Node* n1, const Node* n2) {
        if ((*n1).corePosEdges == (*n2).corePosEdges) {
            if ((*n1).coreTotEdges == (*n2).coreTotEdges) {
                return (*n1).id < (*n2).id;
            }
            return ((*n1).coreTotEdges < (*n2).coreTotEdges);
        }
        return ((*n1).corePosEdges < (*n2).corePosEdges);
    }
    
    bool operator()(const Node &n1, const Node* n2) {
        if (n1.corePosEdges == (*n2).corePosEdges) {
            if (n1.coreTotEdges == (*n2).coreTotEdges) {
                return n1.id < (*n2).id;
            }
            return (n1.coreTotEdges < (*n2).coreTotEdges);
        }
        return (n1.corePosEdges < (*n2).corePosEdges);
    }
            
    bool operator()(const Node* n1, const Node &n2) {
        if ((*n1).corePosEdges == n2.corePosEdges) {
           if ((*n1).coreTotEdges == n2.coreTotEdges) {
                return (*n1).id < n2.id;
            }
            return ((*n1).coreTotEdges < n2.coreTotEdges);
        }
        return ((*n1).corePosEdges < n2.corePosEdges);
    }
};

struct PosCoreComp {
    bool operator()(const Node &n1, const Node &n2) {
        if (n1.positiveCoreNumber == n2.positiveCoreNumber) {
            return n1.id < n2.id;
        }
        return n1.positiveCoreNumber > n2.positiveCoreNumber;
    }
    
    bool operator()(const Node* &n1, const Node* &n2) {
        if ((*n1).positiveCoreNumber == (*n2).positiveCoreNumber) {
            return (*n1).id < (*n2).id;
        }
        return (*n1).positiveCoreNumber > (*n2).positiveCoreNumber;
    }
    
    bool operator()(const Node* n1, const Node* n2) {
        if ((*n1).positiveCoreNumber == (*n2).positiveCoreNumber) {
            return (*n1).id < (*n2).id;
        }
        return (*n1).positiveCoreNumber > (*n2).positiveCoreNumber;
    }
};

unordered_map<int, Node> nodeMap; //Maps IDs to Nodes
CandidateClique maxFound; //The current largest EPC found

/***
* Calculates the Minimum Positive Values
* Returns the gamma value of epsilon/phi
*/
int calculateMPVs(int cap) {
	double x = 1.0;
	double d = 2.0;
	double gamma;
	while (d <= cap) {
		if ((x/d) >= epsilon) {
			if (((d-x)/d) >= phi) {
				gamma = d;
				mpvMap.insert(make_pair((int) d, (int) x));
				d = d + 1.0;
				break;
			} else {
				x = x - 1.0;
				if ((x/d) >= epsilon && ((d-x)/d) >= phi) {
					gamma = d;
					mpvMap.insert(make_pair((int) d, (int) x));
					d = d + 1.0;
					break;
				}
			}
		}
		x = x + 1.0;
		d = d + 1.0;
	}
	
	while (d <= cap) {
		if ((x/d) >= epsilon) {
			if (((d - x)/d) >= phi) {
				mpvMap.insert(make_pair((int) d, (int) x));
			} else {
				mpvMap.insert(make_pair((int) d, 0));
			}
		} else {
			x = x + 1.0; 
			if ((x/d) >= epsilon) {
				if (((d - x)/d) >= phi) {
					mpvMap.insert(make_pair((int) d, (int) x));
				} else {
					mpvMap.insert(make_pair((int) d, 0));
				}
			}
		}
		d = d + 1.0;
	}
	
	return (int) gamma;
}

/***
* Gets the common nodes between searchsPace and the combined list of posN and negN
*/
vector<int> getCommonNodes(vector<int> &searchsPace, vector<int> &posN, vector<int> &negN) {
	vector<int> intersection;
	vector<int>::iterator sIt = searchsPace.begin();
	vector<int>::iterator pIt = posN.begin();
	vector<int>::iterator nIt = negN.begin();
    int sid;
    int pid;
    int nid;
	while (sIt != searchsPace.end() && (pIt != posN.end() || nIt != negN.end())) {
		if (pIt == posN.end()) {
            nid = *nIt;
            sid = *sIt;
			if (nid < sid) {
				nIt++;
			} else if (nid > sid) {
				sIt++;
			} else {
				intersection.push_back(sid);
				sIt++;
				nIt++;
			}
			continue;
		}
		
		if (nIt == negN.end()) {
            pid = *pIt;
            sid = *sIt;
			if (pid < sid) {
				pIt++;
			} else if (pid > sid) {
				sIt++;
			} else {
				intersection.push_back(sid);
				sIt++;
				pIt++;
			}
			continue;
		}
		
        pid = *pIt;
        nid = *nIt;
        sid = *sIt;
		if (pid < nid) {
			if (pid < sid) {
				pIt++;
			} else if (pid > sid) {
				sIt++;
			} else {
				intersection.push_back(sid);
				sIt++;
				pIt++;
			}
		} else {
			if (nid < sid) {
				nIt++;
			} else if (nid > sid) {
				sIt++;
			} else {
				intersection.push_back(sid);
				sIt++;
				nIt++;
			}
		}
	}
	return intersection;
}

/***
* Gets the common nodes between two sorted lists
*/
vector<int> getCommonNodes(vector<int> &list1, vector<int> &list2) {
	vector<int> intersection;
	vector<int>::iterator l1It = list1.begin();
	vector<int>::iterator l2It = list2.begin();
    int id1;
    int id2;
	while (l1It != list1.end() && l2It != list2.end()) {
		id1 = *l1It;
        id2 = *l2It;
		if (id1 < id2) {
			l1It++;
		} else if (id1 > id2) {
			l2It++;
		} else {
			intersection.push_back(*l1It);
			l1It++;
			l2It++;
		}
	}
	return intersection;
}

/***
* Gets the common nodes between two sorted lists (one int, one Node*)
*/
vector<int> getCommonNodes(vector<int> &list1, vector<Node*> &list2) {
	vector<int> intersection;
	vector<int>::iterator l1It = list1.begin();
	vector<Node*>::iterator l2It = list2.begin();
    
    int id1;
    int id2;
	while (l1It != list1.end() && l2It != list2.end()) {
        id1 = *l1It;
        id2 = (*l2It)->id;
		if (id1 < id2) {
			l1It++;
		} else if (id1 > id2) {
			l2It++;
		} else {
			intersection.push_back(*l1It);
			l1It++;
			l2It++;
		}
	}
	return intersection;
}

/***
* Gets the common nodes between two sorted lists improved by considering 
* maxfoundsize
*/
vector<Node*> imprGetCommonNodes(vector<Node*> &list1, vector<Node*> &list2, int &np) {
    vector<Node*> intersection;
    vector<Node*>::iterator l1It = list1.begin();
    vector<Node*>::iterator l2It = list2.begin();
    
    int id1;
    int id2;
    int pcn1;
    int pcn2;
    while (l1It != list1.end() && l2It != list2.end() && 
           (*l1It)->positiveCoreNumber >= maxFoundSize && 
           (*l2It)->positiveCoreNumber >= maxFoundSize) {
        
        id1 = (*l1It)->id;
        id2 = (*l2It)->id;
        pcn1 = (*l1It)->positiveCoreNumber;
        pcn2 = (*l2It)->positiveCoreNumber;
        
        if (pcn1 > pcn2) {
			l1It++;
		} else if (pcn1 < pcn2) {
			l2It++;
		} else {
            if (id1 < id2) {
                l1It++;
            } else if (id1 > id2) {
                l2It++;
            } else if (id1 == id2) {
                intersection.push_back(*l1It);
                np += (*l1It)->negativePotential;
                l1It++;
                l2It++;
            }
			
		}
    }
    return intersection;
}

/***
* Calculates the PCN of a node given the previously found PCN
*/
int getPCNforNode(Node node, int k) {
    unordered_map<int, int>::iterator mpvIt; 
    if (node.coreTotEdges < (k + 1)) {
        return k;
    }
    for (; k <= node.coreTotEdges; k++) {
        mpvIt = mpvMap.find(k + 1);
        if (mpvIt == mpvMap.end()) {
            
            printf("Epsilon and Phi Too Tight\n");
            exit(1);
        }
        if (node.corePosEdges < mpvIt->second || node.coreTotEdges < (k + 1)) {
            return k;
        }
    }
    return node.coreTotEdges;
}

/***
* Decomposes a search space and assigns each node their correct PCN
*/
void posCoreDecomposition(vector<Node*> &space, int gamma) {
	unsigned i = 0;
	while (i < space.size()) {
		space[i]->corePosEdges = space[i]->activePosNeighbours.size();
		space[i]->coreTotEdges = space[i]->activePosNeighbours.size() + space[i]->activeNegNeighbours.size();
        space[i]->active = true;
		i++;
	}
	
    cout << "Original Size: " << space.size() << "\n";
    unordered_map<int, Node>::iterator counterIt;

    NodeComp nc;
    sort(space.begin(), space.end(), nc);
    int k = gamma - 1;
    
    vector<Node*>::iterator spaceIt = space.begin();
    while (!space.empty()) {
        k = getPCNforNode((**spaceIt), k);
        if (k < gamma) {
            (*spaceIt)->positiveCoreNumber = 0;
        } else {
            (*spaceIt)->positiveCoreNumber = k;
        }
        
        
        for (int n : (*spaceIt)->activePosNeighbours) {
            counterIt = nodeMap.find(n);
            vector<Node*>::iterator current = lower_bound(space.begin() + 1, space.end(), counterIt->second, nc);
            if (counterIt->second.active) {
                counterIt->second.corePosEdges--;
                counterIt->second.coreTotEdges--;
            }
            vector<Node*>::iterator newPosition = lower_bound(space.begin() + 1, current, counterIt->second, nc);
            rotate(space.rend() - (current - space.begin()) - 1, space.rend() - (current - space.begin()), space.rend() - (newPosition - space.begin()));
        }
        
        for (int n : (*spaceIt)->activeNegNeighbours) {
            counterIt = nodeMap.find(n);
            vector<Node*>::iterator current = lower_bound(space.begin() + 1, space.end(), counterIt->second, nc);
            if (counterIt->second.active) {
                counterIt->second.coreTotEdges--;
            }
            vector<Node*>::iterator newPosition = lower_bound(space.begin() + 1, current, counterIt->second, nc);
            rotate(space.rend() - (current - space.begin()) - 1, space.rend() - (current - space.begin()), space.rend() - (newPosition - space.begin()));
        }
        (*spaceIt)->active = false;
        spaceIt = space.erase(spaceIt);
    }
    
}

/***
* Checks if a clique satisfies the epsilon and phi requirements
* Baseline version
*/
bool satisfiesRequirements(CandidateClique &clique) {
	if (clique.size() <= 1) {
		return false;
	}
	
	if (clique.getPhi() < phi) {
		return false;
	}
	
	int individualSize = clique.size() - 1;
	for (auto it : clique.nodes) {
		if (it.second.getEpsilon(individualSize) < epsilon) {
			return false;
		}
	}
	
	return true;
}

/***
* Our baseline EPBK algorithm
*/
void baseBK(CandidateClique gRowing, vector<int> searchsPace) {
	if (gRowing.size() > maxFoundSize) {
		if (satisfiesRequirements(gRowing)) {
			maxFoundSize = gRowing.size();
			maxFound = gRowing;
		}
	}
	
	if (searchsPace.empty()) {
		return;
	}
	
	vector<int>::iterator sIt = searchsPace.begin();
	unordered_map<int, Node>::iterator nIt;
	for(; sIt!= searchsPace.end();) {
		nIt = nodeMap.find(*sIt);
		Node nTemp = nIt->second;
		CandidateClique gRowing2 = gRowing;
		gRowing2.adjustCounters(nTemp);
		baseBK(gRowing2, getCommonNodes(searchsPace, nTemp.activePosNeighbours, nTemp.activeNegNeighbours));
		sIt = searchsPace.erase(sIt);
	}
	return;
}

/***
* Improved EPBK function
*/
void imprBK(CandidateClique gRowing, vector<Node*> searchsPace, bool passable) {
    if (passable) {
        if (gRowing.size() > maxFoundSize) {
            maxFoundSize = gRowing.size();
            maxFound = gRowing;
        }
    }
    
    if (searchsPace.empty()) {
        return;
    }

    vector<Node*>::iterator sIt = searchsPace.begin();
    for (; sIt != searchsPace.end();) {

        CandidateClique gRowing2 = gRowing;     
        bool limitBreaker = true;
        bool newPass = gRowing2.adjustAndCheckCounters((**sIt), maxFoundSize, limitBreaker);
        
        if (!limitBreaker) {
            sIt = searchsPace.erase(sIt);
            continue;
        }
        
        int np = 0;
        vector<Node*> searchsPace2 = imprGetCommonNodes(searchsPace, 
                                                        (*sIt)->activeNeighbours, np);
        if (!gRowing2.checkTotalNP(np)) {
            sIt = searchsPace.erase(sIt);
            continue;
        }
        
        imprBK(gRowing2, searchsPace2, newPass);
        
        sIt = searchsPace.erase(sIt);
    }
    return;
}

/***
* Baseline search function (MEPCS)
*/
void searchBaseline(int query) {
	unordered_map<int, Node>::iterator nIt = nodeMap.find(query);
	if (nIt == nodeMap.end()) {
		printf("Invalid Query Node\n");
		exit(3);
	}
	CandidateClique gRowing;
	CandidateNode qNode;
	qNode.id = nIt->second.id;
	qNode.posEdges = 0;
	gRowing.insertNode(qNode);
	gRowing.negEdges = 0;
	gRowing.totEdges = 0;

	vector<int> searchsPace; //Nodes used for search
	vector<int> checksPace;
	for (int n : nIt->second.posNeighbours) {
		searchsPace.push_back(n);
		nIt->second.activePosNeighbours.push_back(n);
	}
	for (int n : nIt->second.negNeighbours) {
		searchsPace.push_back(n);
		nIt->second.activeNegNeighbours.push_back(n);
	}
	sort(searchsPace.begin(), searchsPace.end());
	checksPace = searchsPace;
	checksPace.push_back(nIt->second.id);
	sort(checksPace.begin(), checksPace.end());
	
	unordered_map<int, Node>::iterator uIt;
	for (int n : nIt->second.activePosNeighbours) {
		uIt = nodeMap.find(n);
		uIt->second.activePosNeighbours = getCommonNodes(uIt->second.posNeighbours, checksPace);
		uIt->second.activeNegNeighbours = getCommonNodes(uIt->second.negNeighbours, checksPace);
	}
	
	for (int n : nIt->second.activeNegNeighbours) {
		uIt = nodeMap.find(n);
		uIt->second.activePosNeighbours = getCommonNodes(uIt->second.posNeighbours, checksPace);
		uIt->second.activeNegNeighbours = getCommonNodes(uIt->second.negNeighbours, checksPace);
	}
	
	cout << "Search Space: " << searchsPace.size() << "\n";
	baseBK(gRowing, searchsPace);
	cout << "Max Found Size: " << maxFoundSize << "\n";
	for (auto it : maxFound.nodes) {
		cout << it.second.id << ": " << it.second.posEdges << "\n";
	}
	cout << "Neg Edges: " << maxFound.negEdges << "\n";
	cout << "Tot Edges: " << maxFound.totEdges << "\n";
}

/***
* Improved search function (IMEPCS)
*/
void searchImpr(int q) {
    unordered_map<int, Node>::iterator nIt = nodeMap.find(q);
    int gamma = calculateMPVs(maxDeg);
    cout << "Gamma: " << gamma << "\n";

	vector<Node*> coreSpace; //Nodes used for search
    vector<int> checksPace;    
    coreSpace.push_back(&nIt->second);
    checksPace.push_back(nIt->second.id);
	for (int n : nIt->second.posNeighbours) {
		nIt->second.activePosNeighbours.push_back(n);
        checksPace.push_back(n);
	}
	for (int n : nIt->second.negNeighbours) {
		nIt->second.activeNegNeighbours.push_back(n);
        checksPace.push_back(n);
	}
	
    sort(checksPace.begin(), checksPace.end());
    
	unordered_map<int, Node>::iterator uIt;
	for (int n : nIt->second.activePosNeighbours) {
		uIt = nodeMap.find(n);
		uIt->second.activePosNeighbours = getCommonNodes(uIt->second.posNeighbours, checksPace);
		uIt->second.activeNegNeighbours = getCommonNodes(uIt->second.negNeighbours, checksPace);
        coreSpace.push_back(&uIt->second);
	}
	
	for (int n : nIt->second.activeNegNeighbours) {
		uIt = nodeMap.find(n);
		uIt->second.activePosNeighbours = getCommonNodes(uIt->second.posNeighbours, checksPace);
		uIt->second.activeNegNeighbours = getCommonNodes(uIt->second.negNeighbours, checksPace);
        coreSpace.push_back(&uIt->second);
	}
    
    posCoreDecomposition(coreSpace, gamma);
    
    if (nIt->second.positiveCoreNumber == 0) {
        cout << "No Possible Solution via Positive Core Decomp\n";
        return;
    }
    
    vector<Node*> searchsPace; //Nodes used for search
	checksPace.clear();
    
    //
    vector<int>::iterator killIt = nIt->second.activePosNeighbours.begin();
    for (; killIt != nIt->second.activePosNeighbours.end(); ) {
        uIt = nodeMap.find(*killIt);
        if (uIt->second.positiveCoreNumber == 0) {
            killIt = nIt->second.activePosNeighbours.erase(killIt);
        } else {
            searchsPace.push_back(&uIt->second);
            killIt++;
        }
    }
    
    killIt = nIt->second.activeNegNeighbours.begin();
    for (; killIt != nIt->second.activeNegNeighbours.end(); ) {
        uIt = nodeMap.find(*killIt);
        if (uIt->second.positiveCoreNumber == 0) {
            killIt = nIt->second.activeNegNeighbours.erase(killIt);
        } else {
            searchsPace.push_back(&uIt->second);
            killIt++;
        }
    }

    PosCoreComp pcc;
    sort(searchsPace.begin(), searchsPace.end(), pcc);
    
    //Build Active Neighbour Lists
    unordered_map<int, Node>::iterator tIt;
	for (int n : nIt->second.activePosNeighbours) {
		uIt = nodeMap.find(n);
        killIt = uIt->second.activePosNeighbours.begin();
        for (; killIt != uIt->second.activePosNeighbours.end(); ) {
            tIt = nodeMap.find(*killIt);
            if (tIt->second.positiveCoreNumber == 0) {
                killIt = uIt->second.activePosNeighbours.erase(killIt);
            } else {
                uIt->second.activeNeighbours.push_back(&tIt->second);
                killIt++;
            }
        }
        killIt = uIt->second.activeNegNeighbours.begin();
        for (; killIt != uIt->second.activeNegNeighbours.end(); ) {
            tIt = nodeMap.find(*killIt);
            if (tIt->second.positiveCoreNumber == 0) {
                killIt = uIt->second.activeNegNeighbours.erase(killIt);
            } else {
                uIt->second.activeNeighbours.push_back(&tIt->second);
                killIt++;
            }
        }
        sort(uIt->second.activeNeighbours.begin(), uIt->second.activeNeighbours.end(), pcc);
    }
	
	for (int n : nIt->second.activeNegNeighbours) {
		uIt = nodeMap.find(n);
        
        killIt = uIt->second.activePosNeighbours.begin();
        for (; killIt != uIt->second.activePosNeighbours.end(); ) {
            tIt = nodeMap.find(*killIt);
            if (tIt->second.positiveCoreNumber == 0) {
                killIt = uIt->second.activePosNeighbours.erase(killIt);
            } else {
                uIt->second.activeNeighbours.push_back(&tIt->second);
                killIt++;
            }
        }
        killIt = uIt->second.activeNegNeighbours.begin();
        for (; killIt != uIt->second.activeNegNeighbours.end(); ) {
            tIt = nodeMap.find(*killIt);
            if (tIt->second.positiveCoreNumber == 0) {
                killIt = uIt->second.activeNegNeighbours.erase(killIt);
            } else {
                uIt->second.activeNeighbours.push_back(&tIt->second);
                killIt++;
            }
        }
        sort(uIt->second.activeNeighbours.begin(), uIt->second.activeNeighbours.end(), pcc);
    }
    
    //Calculate Negative Potentials
    for (Node* nPointer : searchsPace) {
        (*nPointer).calculateNP();
    }
    nIt->second.calculateNP();
    
    CandidateClique gRowing;
	CandidateNode qNode;
	qNode.id = nIt->second.id;
	qNode.posEdges = 0;
    qNode.negativeLimit = nIt->second.negativeLimit;
	gRowing.insertNode(qNode);
	gRowing.negEdges = 0;
	gRowing.totEdges = 0;
    gRowing.currentNP = nIt->second.negativePotential;
    
    
	cout << "Post Prune Search Space: " << searchsPace.size() << "\n";
    imprBK(gRowing, searchsPace, false);
    cout << "Max Found Size: " << maxFoundSize << "\n";
    for (auto it : maxFound.nodes) {
		cout << it.second.id << ": " << it.second.posEdges << "\n";
	}
	cout << "Neg Edges: " << maxFound.negEdges << "\n";
	cout << "Tot Edges: " << maxFound.totEdges << "\n";
}

/***
* Initialises the relevent data structures from the Edge File
*
* Edge File Format ((node1)\t(node2)\t(sign))
* sign = 1 if pos and sign = 0 if neg
*/
void initialiseEdges(char* filename) {
	ifstream nameStream(filename);
	if (!nameStream.is_open()) {
		printf("Edges File DNE\n");
		exit(3);
	}
	char line[40];
	char* tempV;
	int v1, v2, sign;
	unordered_map<int, Node>::iterator nIt1, nIt2;
    vector<int> vIDs;

	int posCount = 0;
	int negCount = 0;
	//int sEdgeCount = 0;
	while (nameStream.good()) {
		nameStream.getline(line, 40);
		if (nameStream.eof()) {
			nameStream.close();
			break;
		}
		
		tempV = strtok(line, "\t");
		v1 = stoi(tempV);
		tempV = strtok(NULL, "\t");
		v2 = stoi(tempV);
		tempV = strtok(NULL, "\t");
		sign = stoi(tempV);
		
		nIt1 = nodeMap.find(v1);
		if (nIt1 == nodeMap.end()) {
			Node n;
			n.id = v1;
			nodeMap.insert(make_pair(v1, n));
			nIt1 = nodeMap.find(v1);
			vIDs.push_back(v1);
		}
		
		nIt2 = nodeMap.find(v2);
		if (nIt2 == nodeMap.end()) {
			Node n;
			n.id = v2;
			nodeMap.insert(make_pair(v2, n));
			nIt2 = nodeMap.find(v2);
			vIDs.push_back(v2);
		}
		
		if (sign == 1) {
			nIt1->second.posNeighbours.push_back(v2);
			nIt2->second.posNeighbours.push_back(v1);
			posCount++;
		} else {
			nIt1->second.negNeighbours.push_back(v2);
			nIt2->second.negNeighbours.push_back(v1);
			negCount++;
		}		
		
	}
	
	sort(vIDs.begin(), vIDs.end());
    maxDeg = 0;
	for (int i : vIDs) {
		nIt1 = nodeMap.find(i);
		sort(nIt1->second.posNeighbours.begin(), nIt1->second.posNeighbours.end());
		sort(nIt1->second.negNeighbours.begin(), nIt1->second.negNeighbours.end());
        
        if ((nIt1->second.posNeighbours.size() + nIt1->second.negNeighbours.size()) > maxDeg) {
            maxDeg = nIt1->second.posNeighbours.size() + nIt1->second.negNeighbours.size();
        }
	}
	cout << "Nodes: " << vIDs.size() << "\n";
	cout << "Total Edges: " << posCount + negCount << "\n";
	cout << "Positive Edges: " << posCount << "\n";
	cout << "Negative Edges: " << negCount << "\n";
}

/***
* Initialises the list of queries to conduct
*/
void initialiseQueryList(char* filename, vector<int> &queryList) {
	ifstream nameStream(filename);
	if (!nameStream.is_open()) {
		printf("Edges File DNE\n");
		exit(3);
	}
	char line[40];
	char* tempV;
	int v;
	while (nameStream.good()) {
		nameStream.getline(line, 40);
		if (nameStream.eof()) {
			nameStream.close();
			break;
		}
		
		tempV = strtok(line, "\t");
		v = stoi(tempV);
		queryList.push_back(v);
	}
	
	cout << "# of Queries: " << queryList.size() << "\n";
}

/***
* argv[0] = epcli
* argv[1] = function
		1: Single Search Baseline (MEPCS)
        2: Single Search Improved (IMEPCS)
* argv[2] = EDGELIST
* argv[3] = Query Node/QUERYFILE
* argv[4] = Epsilon Value (Range [0, 1])
* argv[5] = Phi Value (Range [0, 1])
* argv[6] = OUTFILE
*/
int main(int argc, char *argv[]) {
	if (argc < 4 || argc > 7) {
		printf("Insufficient Args\n");
		exit(1);
	}
	epsilon = stod(argv[4]);
	phi = stod(argv[5]);	
    maxFoundSize = 0;
	int function = stoi(argv[1]);
	cout << "Loading Edges\n";
	initialiseEdges(argv[2]);
	auto start = chrono::high_resolution_clock().now();
	
	//Functions
	if (function == 1) {
		searchBaseline(stoi(argv[3]));
	} else if (function == 2) {
        searchImpr(stoi(argv[3]));  
    } else {
		printf("Invalid Function Call\n");
		exit(1);
	}
	
	auto stop = chrono::high_resolution_clock().	now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start).count();
	cout << "Total Execution Time: " << duration << "\n";
    if (argc == 7) {
        string writeString = "Total Runtime: ";
        writeString += to_string(duration);
        writeString += "\nSize: ";
        writeString += to_string(maxFoundSize);
        writeString += "\n";
        writeToFile(writeString, argv[6]);
    }

}

