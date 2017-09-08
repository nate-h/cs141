#ifndef HEADERFILE_H
#define HEADERFILE_H


//includes
#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <fstream>
#include <sstream>



int myMax(int val1, int val2) {   // renamed because scoping
	if (val1 >= val2)
		return val1;
	return val2;
}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////    GRAPH STRUCTURES   //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct Edge;

struct Node {
	int value;
        int refPos;
	std::vector<Edge*>* neighbors;

  Node(int v,int r) :
    value(v), refPos(r) {
	  neighbors = new std::vector<Edge*>();
	}
	void addEdge(Edge * e) {
		neighbors->push_back(e);
	}
	~Node() {
	}
};

struct Edge {
	Node * end;
	int cost;
	char how;
	Edge(Node * en, int c, char h) :
			end(en), cost(c), how(h) {
	}
	~Edge() {
	}
};

struct traceBackNodes
{
  int ct;
  int val;
  int posY;
  Node * current;
  std::string r;
  std::string q;
  traceBackNodes(int c, int v, int y, Node * start, std::string a, std::string b) :
    ct(c), val(v), posY(y),
    current(start),
    r(a), q(b)
  {}
};

struct listVal {
	Node* score;
	listVal* next;
	listVal* prev;
	listVal(Node* addMe) {
		score = addMe;
		next = NULL;
		prev = NULL;
	}
};


struct Graph {
	std::vector <Node*>* graph;
	Graph() {
		graph = new std::vector<Node*>();
	}
	~Graph() {
	  delete graph;
	}
  void addNode(int val, int r) {
    Node * temp = new Node(val,r);
    graph->push_back(temp);
	}
	void addEdge(Node *n1, Node *n2, int cst, char hw) {
		Edge *e = new Edge(n2, cst, hw);
		n1->addEdge(e);
	}

};

////////////////////////////////////////////////////////////////////////////////
/////////////////////    Score Calculating Function   //////////////////////////
////////////////////////////////////////////////////////////////////////////////

//switch cases are faster than nested if else statements
//jump tables/ look up tables
int scores(char ref, char qry) {

	switch (ref) {
	case 'A':
		switch (qry) {
		case 'A':
			return 100;
			break;
		case 'G':
			return -10;
			break;
		case 'T':
			return -15;
			break;
		case 'C':
			return -10;
			break;
		default:
			std::cout << "THIS SHOULDNT PRINT!" << std::endl;
			break;
		}
		break;
	case 'G':
		switch (qry) {
		case 'A':
			return -10;
			break;
		case 'G':
			return 100;
			break;
		case 'T':
			return -10;
			break;
		case 'C':
			return -15;
			break;
		default:
			std::cout << "THIS SHOULDNT PRINT!" << std::endl;
			break;
		}
		break;
	case 'T':
		switch (qry) {
		case 'A':
			return -15;
			break;
		case 'G':
			return -10;
			break;
		case 'T':
			return 100;
			break;
		case 'C':
			return -10;
			break;
		default:
			std::cout << "THIS SHOULDNT PRINT!" << std::endl;
			break;
		}
		break;
	case 'C':
		switch (qry) {
		case 'A':
			return -10;
			break;
		case 'G':
			return -15;
			break;
		case 'T':
			return -10;
			break;
		case 'C':
			return 100;
			break;
		default:
			std::cout << "THIS SHOULDNT PRINT!" << std::endl;
			break;
		}
		break;
	default:
		std::cout << "THIS SHOULDNT PRINT!" << std::endl;
		break;
	}
	return 0;

}

#endif
