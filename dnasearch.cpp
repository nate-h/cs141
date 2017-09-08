////////////////////////////////////////////////////////////////////////////////
//  Project 3: DNA Search
//  Team Members:
//  Stephen    sbola004@ucr.edu.com
//  Mike       mrome007@ucr.edu.com
//  Jon        jdean007@ucr.edu.com
//  Nate         nhapeman@gmail.com
//  Aaron       aaronmg83@gmail.com
////////////////////////////////////////////////////////////////////////////////

//includes
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include "definitions.h"
#include <utility>
#include <stdlib.h>

//forward includes
int scores(char ref, char qry);
bool readInputFiles(char* dataBaseName, char* queryName);

//globals
int key = 0;
char * query;
char* dataBase;
char * finalQuery;
char* finalDataBase;
int dataBaseLength;
int queryLength;
int num = 0;
int nodeToDel = 0;
int nS;
////////////////////////////////////////////////////////////////////////////////
/////////////////////////////   FIND ALIGNMENT   ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void maxAlign(Graph &g, std::vector<Node*> &lastRow, std::vector<Node*> topRow,
		std::vector<Node*> fcol, char * align) {

	std::vector<Node*> check(fcol.size());
	int bottom = queryLength - 1;
	lastRow[0] = fcol[bottom];
	for (unsigned i = 1; i < topRow.size(); i++) {
		check[0] = topRow[i];
		for (unsigned j = 1; j < fcol.size(); j++) {
			int gp_vertical = 0;
			int gp_diagonal = 0;
			int gp_horizontal = 0;
			gp_diagonal = scores(dataBase[i], query[j]);

			bool hasVert = false;
			bool hasHori = false;
			for (unsigned v = 0; v < check[j - 1]->neighbors->size(); v++) {
				if (check[j - 1]->neighbors->at(v)->how == 'v')
					hasVert = true;
			}

			for (unsigned v = 0; v < fcol[j]->neighbors->size(); v++) {
				if (fcol[j]->neighbors->at(v)->how == 'h')
					hasHori = true;
			}

			if (hasVert == false)
				gp_vertical = -20;
			else
				gp_vertical = -5;
			if (hasHori == false)
				gp_horizontal = -20;
			else
				gp_horizontal = -5;

			int ver = check[j - 1]->value + gp_vertical;
			int dia = fcol[j - 1]->value + gp_diagonal;
			int hor = fcol[j]->value + gp_horizontal;

			int mx = myMax(dia, myMax(ver, hor));
			g.addNode(mx, i);
			check[j] = g.graph->at(key);

			if (dia > ver && dia > hor)
				g.addEdge(check[j], fcol[j - 1], gp_diagonal, 'd');
			else if (ver > dia && ver > hor)
				g.addEdge(check[j], check[j - 1], gp_vertical, 'v');
			else if (hor > dia && hor > ver)
				g.addEdge(check[j], fcol[j], gp_horizontal, 'h');
			else if (dia > ver && dia == hor) {
				if (hasVert)
					g.addEdge(check[j], fcol[j - 1], gp_diagonal, 'd');
				g.addEdge(check[j], fcol[j], gp_horizontal, 'h');
			} else if (dia == ver && dia > hor) {
				if (hasHori)
					g.addEdge(check[j], fcol[j - 1], gp_diagonal, 'd');
				g.addEdge(check[j], check[j - 1], gp_vertical, 'v');
			} else if (hor > dia && hor == ver) {
				g.addEdge(check[j], check[j - 1], gp_vertical, 'v');
				g.addEdge(check[j], fcol[j], gp_horizontal, 'h');
			} else if (dia == ver && dia == hor) {
				g.addEdge(check[j], fcol[j - 1], gp_diagonal, 'd');
				g.addEdge(check[j], fcol[j], gp_horizontal, 'h');
				g.addEdge(check[j], check[j - 1], gp_vertical, 'v');
			}
			key++;
		}

		for (unsigned k = 0; k < fcol.size() - 1; k++) {
			if (align[0] == '0') {
				if (k != fcol.size() - 2) {
					for (unsigned l = 0;
							l < g.graph->at(nodeToDel)->neighbors->size(); l++)
						delete g.graph->at(nodeToDel)->neighbors->at(l);
					delete g.graph->at(nodeToDel)->neighbors;
					delete g.graph->at(nodeToDel);
				}
				nodeToDel++;
			}
			fcol[k] = check[k];
		}
		fcol[fcol.size() - 1] = check[fcol.size() - 1];
		lastRow[i] = check[bottom];
	}
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////   FINDING A PATH   ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void thePaths(Node * start, const char * ref, const char * que) {
	int quePos = queryLength - 1;
	std::stack<traceBackNodes*> toVisit;
	Node * first = start;
	toVisit.push(new traceBackNodes(0, start->value, quePos, first, "", ""));
	while (!toVisit.empty()) {
		int rx = toVisit.top()->current->refPos;
		int qx = toVisit.top()->posY;
		Node * curr = toVisit.top()->current;
		std::string currR = toVisit.top()->r;
		std::string currQ = toVisit.top()->q;
		int va = toVisit.top()->val;
		int cos = toVisit.top()->ct;
		if (curr->neighbors->size() == 0) {
		  double sc = (start->value - va)/100.0;
			std::cout << "score: " << sc << std::endl;
			std::cout << currR << std::endl;
			std::cout << currQ << std::endl;
			std::cout << std::endl;
			num++;
			if (num > nS-2) //change number here to change the number of alignments
				break;    //to be printed
		}
		Edge * temp = NULL;
		//std::cout << va << std::endl;
		delete toVisit.top();
		toVisit.pop();
		int c = 0;
		for (unsigned i = 0; i < curr->neighbors->size(); i++) {
			temp = curr->neighbors->at(i);
			Node * nei = temp->end;
			if (temp->how == 'd') {
				int y = qx - 1;
				std::string a = ref[rx] + currR;
				std::string b = que[qx] + currQ;
				toVisit.push(
						new traceBackNodes(temp->cost, va - temp->cost, y, nei,
								a, b));
			} else if (temp->how == 'h') {
				int y = qx;
				std::string a = ref[rx] + currR;
				std::string b = "_" + currQ;

				if (cos == -20 || cos == -5)
					c = -5;
				else
					c = -20;

				toVisit.push(new traceBackNodes(c, va - c, y, nei, a, b));
			} else if (temp->how == 'v') {
				int y = qx - 1;
				std::string a = "_" + currR;
				std::string b = que[qx] + currQ;
				if (cos == -20 || cos == -5)
					c = -5;
				else
					c = -20;

				toVisit.push(new traceBackNodes(c, va - c, y, nei, a, b));
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////////   FINDING TOP VALUES   //////////////////////////////
////////////////////////////////////////////////////////////////////////////////

std::vector<Node*> topValues(std::vector<Node*> topRow, int howMany) {

	if (topRow.size() == 0) {
		std::cout << "Error: You passed a vector of size zero!" << std::endl;
	}

	//initializing place holders
	listVal* head = new listVal(topRow[0]);
	listVal* tail = head;

	//iterators
	listVal* iter = head;
	bool used = false;

	//------------------------------------------------ add first "howMany" nodes
	for (int i = 1; i < howMany; i++) {

		int val = topRow[i]->value;
		listVal* addMe = new listVal(topRow[i]);
		iter = head;

		//replace head if val > max
		if (val > head->score->value) {
			head = addMe;
			head->next = iter;
			iter->prev = head;
		} else {

			used = false;
			if (iter->next != NULL)
				iter = iter->next;

			while (iter->next != NULL) {

				//before
				if (val >= iter->score->value) {
					addMe->prev = iter->prev;
					addMe->next = iter;
					iter->prev->next = addMe;
					iter->prev = addMe;
					used = true;
					break;
				}
				iter = iter->next;
			}

			if (used == false) {
				if (val > tail->score->value) {
					addMe->prev = tail->prev;
					tail->prev->next = addMe;
					tail->prev = addMe;
					addMe->next = tail;

				} else {
					iter = tail;
					addMe->prev = tail;
					tail = addMe;
					iter->next = addMe;
				}

			}
		}

	}

	//--------------------------------------------- adding the rest of the nodes

	for (unsigned i = howMany; i < topRow.size(); i++) {

		int val = topRow[i]->value;

		if (val > tail->score->value) {

			iter = head;

			//replace head if val > max
			if (val > head->score->value) {
				listVal* addMe = new listVal(topRow[i]);
				head = addMe;
				head->next = iter;
				iter->prev = head;

				//deleting off excess
				iter = tail;
				tail->prev->next = NULL;
				tail = tail->prev;
				//delete iter->score;
				delete iter;
			} else {
				listVal* addMe = new listVal(topRow[i]);
				//used = false;
				if (iter->next != NULL)
					iter = iter->next;

				while (iter->next != NULL) {

					//before
					if (val >= iter->score->value) {
						addMe->prev = iter->prev;
						addMe->next = iter;
						iter->prev->next = addMe;
						iter->prev = addMe;
						//deleting off excess
						iter = tail;
						tail->prev->next = NULL;
						tail = tail->prev;
						//delete iter->score;
						delete iter;
						break;
					}
					iter = iter->next;
				}
			}

		}
	}
	//return vector of largest nodes
	std::vector<Node*> returnMe;
	iter = head;
	while (iter != NULL) {
		returnMe.push_back(iter->score);
		iter = iter->next;
	}

	return returnMe;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////   MAIN FUNCTION   /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

	//used to store command line arguments
	char*dataBaseName;
	char*queryName;
	char*numSequences;
	char*alignment;
	
	//parsing through input arguments
	for (int i = 0; i < argc; ++i) {
		char*temp = argv[i];
		if (strcmp(temp, "d") == 0)
			dataBaseName = argv[i + 2];
		else if (strcmp(temp, "q") == 0)
			queryName = argv[i + 2];
		else if (strcmp(temp, "n") == 0)
			numSequences = argv[i + 2];
		else if (strcmp(temp, "a") == 0)
			alignment = argv[i + 2];
	}

	//testing validity
	std::cout << "Here are the input arguments" << std::endl;
	std::cout << "DataBase: " << dataBaseName << std::endl;
	std::cout << "Query: " << queryName << std::endl;
	std::cout << "Number Sequences: " << numSequences << std::endl;
	std::cout << "Alignment: " << alignment << std::endl;

	readInputFiles(dataBaseName, queryName);
	nS = atoi(numSequences);
	//dataBase = referPlusSpace;
	//query = queryPlusSpace;
	/////////////////////////////////////////////////////
	if(nS > dataBaseLength)
	  {
	    std::cout << "Not Enough Sequences" << std::endl;
	    return 0;
	  }
	std::cout << "DataBase Length: " << dataBaseLength << std::endl;
	std::cout << "Query Length: " << queryLength << std::endl;
	//top row
	std::vector<Node*> topRow(dataBaseLength);
	//first col
	std::vector<Node*> firstCol(queryLength);

	Graph g;

	for (int i = 0; i < dataBaseLength; i++) {
		g.addNode(0, i);
		topRow[i] = g.graph->at(key);
		key++;
	}

	nodeToDel = key;
	firstCol[0] = topRow[0];
	for (int i = 1; i < queryLength; i++) {
		if (firstCol[i - 1]->neighbors->size() == 0) {
			int v = -20 + firstCol[i - 1]->value;
			g.addNode(v, 0);
			firstCol[i] = g.graph->at(key);
			g.addEdge(firstCol[i], firstCol[i - 1], -20, 'v');
		} else {
			int v = -5 + firstCol[i - 1]->value;
			g.addNode(v, 0);
			firstCol[i] = g.graph->at(key);
			g.addEdge(firstCol[i], firstCol[i - 1], -5, 'v');
		}
		key++;
	}

	std::vector<Node*> lastRow(dataBaseLength);
	maxAlign(g, lastRow, topRow, firstCol, alignment);
	std::cout << std::endl;
	std::cout << "============================================" << std::endl;
	//the number depends on number sequences
	std::vector<Node*> sortedVec = topValues(lastRow, nS);

	if (alignment[0] == '0') {
		//the number depends on number sequences;
	  double sc = 0;
		for (int i = 0; i < nS; i++)
		  {
		    sc = sortedVec[i]->value/100.0;
		    std::cout << sc << " ";
		  }
		std::cout << std::endl;
	} else {
		//depends on number sequences
		for (int i = 0; i < nS && num < nS; i++) {
			std::cout << "============================================"
					<< std::endl;
			std::cout << "The alignments are: " << std::endl;
			thePaths(sortedVec[i], dataBase, query);
		}
	}

	////////////////////////////////////////////////////////////////

	// deleting dynamically allocated memory
	delete[] dataBase;
	delete[] query;
	if (alignment[0] == '1') {
		for (int i = g.graph->size() - 1; i >= 0; i--) {
			for (unsigned j = 0; j < g.graph->at(i)->neighbors->size(); j++) {
				delete g.graph->at(i)->neighbors->at(j);
			}
			delete g.graph->at(i)->neighbors;
			delete g.graph->at(i);
		}
	}

	return 0;

}

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////   File Reading   /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool readInputFiles(char* dataBaseName, char* queryName) {

	//--------------------------------------------------------- loading dataBase

	std::ifstream dataFile(dataBaseName, std::ifstream::binary);

	//transferring data to database pointer
	if (dataFile) {

		dataFile.seekg(0, dataFile.end);
		dataBaseLength = dataFile.tellg();
		dataFile.seekg(0, dataFile.beg);
		dataBase = new char[dataBaseLength + 1];
		finalDataBase = new char[dataBaseLength + 1];
		dataFile.read(dataBase, dataBaseLength);
		dataFile.close();
		dataBase[dataBaseLength] = '\0';
		int skipHowMuch=1;
		if(dataBase[0]=='>')
		{
			for(int i =0; i< dataBaseLength ; ++i)
			{
				if(dataBase[i] == '\n')
					break;
				skipHowMuch++;

			}

		}

		int count=1;
		finalDataBase[0]=' ';
		for(int i =skipHowMuch; i < dataBaseLength; ++i){
			if(dataBase[i] == 'T' || dataBase[i] == 'A' || dataBase[i] == 'G' || dataBase[i] == 'C')
			{finalDataBase[count]=dataBase[i]; ++count;}
		}

		finalDataBase[count]='\0';
		dataBaseLength=count;
		delete []dataBase;
		dataBase=finalDataBase;
		std::cout<<"DataBase"<< finalDataBase<<std::endl;

	} else {
		std::cout << "Error: Where is dataBase" << std::endl;
		return false;
	}

	//------------------------------------------------------------ loading query

	std::ifstream queryFile(queryName, std::ifstream::binary);

	//transferring data to database pointer
	if (queryFile) {

		queryFile.seekg(0, queryFile.end);
		queryLength = queryFile.tellg();
		queryFile.seekg(0, queryFile.beg);
		query = new char[queryLength + 1];
		finalQuery = new char[queryLength + 1];
		queryFile.read(query, queryLength);
		queryFile.close();
		query[queryLength] = '\0';

		int count=1;
		finalQuery[0]=' ';
		for(int i =0; i < queryLength; ++i){
			if(query[i] == 'T' || query[i] == 'A' || query[i] == 'G' || query[i] == 'C')
			{finalQuery[count]=query[i]; ++count;}
		}

		finalQuery[count]='\0';
		queryLength=count;
		delete []query;
		query=finalQuery;
		std::cout<<"query"<< finalQuery<<std::endl;

	} else {
		std::cout << "Error: Where is query" << std::endl;
		return false;
	}

	return true;
}

