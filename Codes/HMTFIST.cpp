//Assumption no missing Value in this test data
#define Dim 3 //input data dimension
#define cNum 2
#define _USE_MATH_DEFINES
#define LOGZERO -INFINITY
#define _CRT_SECURE_NO_WARNINGS
#define MESSAGELOW -INFINITY
#define MESSAGEHIGH 10000
#define PIXELLIMT 3000
#define MAXCOST -1000.0

#define EB 2  //1 use max cost for identifying effective branches 0 use chain length
#define BOUNDARY_NODES_OBSERVED 0   // 0 does not exclude any branches // 1 consider pits and tree and unobserved and exclude those branches // 2 excludes branches if the boundary of flood and dry is overlapping with pits layer
#define EFFECTIVE_BRANCH_VIZ 1
#define DEBUG_OUTPUT 1
#include<iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include<fstream>
#include<algorithm>
#include<numeric>
#include<vector>
#include<string>
#include<chrono>
#include<ctime>
#include<cmath>
#include<limits>
#include<cstdio>
#include<queue>
#include <stack> 
#include <list>
#include<unordered_set> 
#include <iomanip>
#include <sstream>
#include <map>
#include "GeotiffRead.cpp"
#include "GeotiffWrite.cpp"
#include "Tree.cpp"
#include "DataTypes.h"
using namespace std;


class cFlood {
private:
	//struct sComponent unionComponent;
	struct sParameter parameter;
	struct sData data;
	vector<struct sData> subtreeData;
	struct sTree tree;
	struct sInference infer;
	vector<int>testIndex;
	vector<int>testLabel;
	vector<int>mappredictions;
	//ofstream timeLogger;
	std::string CTInputLocation;
	std::string CTSourceDirection;
	std::string CTProbability;
	std::string CTBank;
	std::string CTCost;
	std::string CTPits;
	std::string CTTree;
	std::string CTRoughness;
	std::string CTFel;
	std::string CTPara;
	std::string CTStream;
	std::string CTOutputLocation;
	std::string CTPrediction;

	//tree construction
	struct subset* subsets;

	//new
	vector<double>elnPzn_xn;

public:
	void input(int argc, char* argv[]);


	void UpdateTransProb(); //Update P(y|z), P(zn|zpn), P(zn|zpn=empty)

	//inference 
	void inference();
	void interpolate();
	void output();
	void prediction();

	//helper functions
	void removeLink(vector<int>& v, int removeID);
	void displayTree(int TreeID);
	void updateMapPrediction_left();
	void updateMapPrediction_right();
	vector<int> getBFSOrder(int root, vector<int>& bfsVisited, int bank);
	//struct conMatrix getConfusionMatrix();

	//Test Modules
	void sanityChecker();
	void getOriginIdBanks();
	void getOrgIds();
	void getIds();
	void getOriginIdLeftBanks();
	void getOriginIdRightBanks();
	void getRegionNodeCount();
	void getLeftRegionNodeCount();
	void getRightRegionNodeCount();
	void getOriginIdBanks_effectiveBranches();
};



void getCofactor(double mat[Dim][Dim], double temp[Dim][Dim], int p, int q, int n) {
	int i = 0, j = 0;
	// Looping for each element of the matrix
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
			//  Copying into temporary matrix only those element
			//  which are not in given row and column
			if (row != p && col != q) {
				temp[i][j++] = mat[row][col];

				// Row is filled, so increase row index and
				// reset col index
				if (j == n - 1) {
					j = 0;
					i++;
				}
			}
		}
	}
}

//dynamic memory allocation,dimensional two dimension array
/* Recursive function for finding determinant of matrix.
n is current dimension of mat[][]. */
double determinant(double mat[Dim][Dim], int n) {
	double D = 0; // Initialize result

				  //  Base case : if matrix contains single element
	if (n == 1)
		return mat[0][0];

	double temp[Dim][Dim]; // To store cofactors
	int sign = 1;  // To store sign multiplier

				   // Iterate for each element of first row
	for (int f = 0; f < n; f++) {
		// Getting Cofactor of mat[0][f]
		getCofactor(mat, temp, 0, f, n);
		D += sign * mat[0][f] * determinant(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}
	return D;
}

// Function to get adjoint of A[Dim][Dim] in adj[Dim][Dim].
void adjoint(double A[Dim][Dim], double adj[Dim][Dim]) {
	if (Dim == 1) {
		adj[0][0] = 1;
		return;
	}

	// temp is used to store cofactors of A[][]
	int sign = 1;
	double temp[Dim][Dim];

	for (int i = 0; i < Dim; i++) {
		for (int j = 0; j < Dim; j++) {
			// Get cofactor of A[i][j]
			getCofactor(A, temp, i, j, Dim);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j][i] = (sign) * (determinant(temp, Dim - 1));
		}
	}
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(double A[Dim][Dim], double inverse[Dim][Dim]) {
	// Find determinant of A[][]

	if (Dim == 1) {
		inverse[0][0] = 1.0 / A[0][0];
		return true;
	}

	double det = determinant(A, Dim);
	if (det == 0) {
		std::cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint
	double adj[Dim][Dim];
	adjoint(A, adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i = 0; i < Dim; i++)
		for (int j = 0; j < Dim; j++)
			inverse[i][j] = adj[i][j] / double(det);
	return true;
}

// extended ln functions
double eexp(double x) {
	if (x == LOGZERO) {
		return 0;
	}
	else {
		return exp(x);
	}
}

double eln(double x) {
	if (x == 0) {
		return LOGZERO;
	}
	else if (x > 0) {
		return log(x);
	}
	else {
		std::cout << "Negative input error " << x << endl;
		exit(0);
	}
}

double elnsum(double x, double y) {
	if (x == LOGZERO) {
		return y;
	}
	else if (y == LOGZERO) {
		return x;
	}
	else if (x > y) {
		return x + eln(1 + eexp(y - x));
	}
	else {
		return y + eln(1 + eexp(x - y));
	}
}

double elnproduct(double x, double y) {
	if (x == LOGZERO || y == LOGZERO) {
		return LOGZERO;
	}
	else {
		return x + y;
	}
}
int dirExists(const char* const path)
{
	struct stat info;

	int statRC = stat(path, &info);
	if (statRC != 0)
	{
		if (errno == ENOENT) { return 0; } // something along the path does not exist
		if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
		return -1;
	}

	return (info.st_mode & S_IFDIR) ? 1 : 0;
}
void cFlood::UpdateTransProb() {
	if (cNum != 2) {
		std::cout << "cannot handle more than two classes now!" << endl;
		std::exit(1);
	}

	double eln(double);
	parameter.elnPz[0] = eln(1 - eexp(parameter.Pi));
	parameter.elnPz[1] = parameter.Pi;
	parameter.elnPz_zpn[0][0] = eln(1);
	parameter.elnPz_zpn[0][1] = parameter.Epsilon;
	parameter.elnPz_zpn[1][0] = eln(0);
	parameter.elnPz_zpn[1][1] = eln(1 - eexp(parameter.Epsilon));
	if (eexp(parameter.Epsilon) < 0 || eexp(parameter.Epsilon) > 1) {
		std::cout << "Epsilon Error: " << eexp(parameter.Epsilon) << endl;
	}
	if (eexp(parameter.Pi) < 0 || eexp(parameter.Pi) > 1) {
		std::cout << "Pi Error: " << eexp(parameter.Pi) << endl;
	}
	if (eexp(parameter.elnPz_zpn[0][1]) + eexp(parameter.elnPz_zpn[1][1]) != 1) {
		std::cout << "Error computing parameter.elnPz_zpn " << endl;
	}
	if (eexp(parameter.elnPz[0]) + eexp(parameter.elnPz[1]) != 1) {
		std::cout << "Error computing parameter.elnPz " << endl;
	}
}

vector<int> cFlood::getBFSOrder(int root, vector<int>& bfsVisited, int bank) {
	//vector<int> bfsVisited;
	vector<int> bfs;
	queue<int> que;
	que.push(root);

	while (!que.empty()) {
		int currentNode = que.front();
		bfs.push_back(currentNode);
		bfsVisited[currentNode] = 1;
		que.pop();
		for (int i = 0; i < data.allNodes[currentNode]->childrenID.size(); i++) {
			int child = data.allNodes[currentNode]->childrenID[i];
			if (!bfsVisited[child]) {
				if (data.allNodes[child]->bank == 0 || data.allNodes[child]->bank == bank) {
					que.push(child);
				}

			}
		}
		for (int i = 0; i < data.allNodes[currentNode]->parentsID.size(); i++) {
			int parent = data.allNodes[currentNode]->parentsID[i];
			if (!bfsVisited[parent]) {
				if (data.allNodes[parent]->bank == 0 || data.allNodes[parent]->bank == bank) {
					que.push(parent);
				}
				//que.push(parent);
			}
		}
	}
	return bfs;
}

void cFlood::input(int argc, char* argv[]) {

	GDALAllRegister();
	clock_t start_s = clock();
	if (argc > 1) {
		ifstream config(argv[1]);
		string line;
		getline(config, line);
		CTInputLocation = line;  //Input file location 
		getline(config, line);
		CTSourceDirection = line;           //Elevation data file name
		getline(config, line);
		CTProbability = line;
		getline(config, line);
		CTBank = line;
		getline(config, line);
		CTCost = line;
		getline(config, line);
		CTPits = line;
		getline(config, line);
		CTTree = line;
		getline(config, line);
		CTRoughness = line;
		getline(config, line);
		CTFel = line;
		getline(config, line);
		CTPara = line;          //parameter data file name
		getline(config, line);
		CTStream = line;
		getline(config, line);
		CTOutputLocation = line; //oputput location to store the output of HMCT
		getline(config, line);
		CTPrediction = line;    //file name for output prediction data

	}
	else {
		std::cout << "Missing Configuration File!";
	}
	struct stat info;

	//if (stat((CTOutputLocation+"testr\\").c_str(), &info) != 0)
	//	printf("cannot access %s\n", (CTOutputLocation + "testr\\").c_str());
	//else if (info.st_mode & S_IFDIR)  // S_ISDIR() doesn't exist on my windows 
	//	printf("%s is a directory\n", (CTOutputLocation + "testr\\").c_str());
	//else
	//	printf("%s is no directory\n", (CTOutputLocation + "testr\\").c_str());
	int status = dirExists(CTInputLocation.c_str());
	if (status <= 0) {
		cout << "Error: input directory does not exist.." << endl;
		exit(0);
	}
	status = dirExists(CTOutputLocation.c_str());
	if (status <=0) {
		cout << "Error: output directory does not exist.." << endl;
		exit(0);
	}
	status = dirExists((CTOutputLocation + "ProfileTables").c_str());
	if (status <= 0) {
		status = _mkdir((CTOutputLocation + "ProfileTables").c_str());
		if (status != 0) {
			cout << "Error: could not create ProfileTables folder.." << endl;
			exit(0);
		}
		//exit(0);
	}

	//reading text file
	ifstream parameterFile(CTInputLocation + CTPara);
	if (!parameterFile) {
		std::cout << "Failed to open parameter!" << endl;
		exit(0);
	}
	parameterFile >> parameter.reachId;
	parameterFile >> parameter.Epsilon;
	parameterFile >> parameter.Pi;
	parameterFile >> parameter.cutoff_percentage;
	parameterFile >> parameter.minCost;

	//parameter.reachId = CTPara.substr(9, 2);
	parameterFile.close();

	//reading stream csv file
	ifstream streamFile(CTInputLocation + CTStream);
	string line;
	getline(streamFile, line); //skipping first line in csv file;
	while (getline(streamFile, line)) {
		//std::cout << line<<endl;
		vector<float> result;
		stringstream s_stream(line);
		while (s_stream.good()) {      //result 0-currentid, 1-sourceDirection 2-bank, 3- observation indicator 4-cost, 5-probability 
			string substr;
			getline(s_stream, substr, ',');
			result.push_back(stof(substr));
		}
		data.reach_ids.push_back((int)result[0]);
		data.distance_covered.push_back((int)result[1]);
		data.flow_accumulation.push_back(result[2]);
	}
	// create Geotiff object for rasters
	//GeotiffRead sourceDirTiff("C:\\Users\\asainju\\Arpan\\IRONFIST\\Data\\omaha_ne_13_fist\\fist6\\19_srcdir.tif");

	GeotiffRead sourceDirTiff((CTInputLocation + CTSourceDirection).c_str());
	GeotiffRead floodProbTiff((CTInputLocation+ CTProbability).c_str());
	GeotiffRead BankTiff((CTInputLocation + CTBank).c_str());
	GeotiffRead CostTiff((CTInputLocation + CTCost).c_str());
	GeotiffRead PitsTiff((CTInputLocation + CTPits).c_str());
	GeotiffRead TreeTiff((CTInputLocation + CTTree).c_str());
	GeotiffRead RoughnessTiff((CTInputLocation + CTRoughness).c_str());
	GeotiffRead FelTiff((CTInputLocation + CTFel).c_str());



	// The pointer to the raster band data of the source direction tiff  
	float** sourceDirData = (float**)sourceDirTiff.GetRasterBand(1); 
	float** floodProbData = floodProbTiff.GetRasterBand(1);
	float** bankData = BankTiff.GetRasterBand(1);
	float** costData = CostTiff.GetRasterBand(1);
	float** pitsData = PitsTiff.GetRasterBand(1);
	float** treeData = TreeTiff.GetRasterBand(1);
	float** roughnessData = RoughnessTiff.GetRasterBand(1);
	float** felData = FelTiff.GetRasterBand(1);

	// Get the array dimensions
	int* dims = sourceDirTiff.GetDimensions();
	// Get the nodata value
	//float noDataValue = (double)sourceDirTiff.GetNoDataValue();

	parameter.ROW = dims[0];
	parameter.COLUMN = dims[1];
	parameter.allPixelSize = parameter.ROW * parameter.COLUMN;
	parameter.orgPixelSize = parameter.allPixelSize;

	//Note end

	if (parameter.Epsilon > 1 || parameter.Pi > 1) {
		std::cout << "wrong parameter" << endl;
	}

	std::cout << "Input parameters:" << endl << "Epsilon: " << parameter.Epsilon << " Pi: " << parameter.Pi << endl;

	// Fill in the data values in the tree
	int nodeIdCounter = 0;
	std::map<int, int> nodeIndexMap;
	int NCOLS = dims[1];
	int NROWS = dims[0];
	for (int row = 0; row < NROWS; row++)
	{
		for (int col = 0; col < NCOLS; col++)
		{
			int sourceDir = sourceDirData[row][col];
			int floodProb = floodProbData[row][col];
			int bank = bankData[row][col];
			int cost = costData[row][col];
			int pits = (pitsData[row][col]==1)? 1: 0;
			int tree = (treeData[row][col]>0)? 1: 0;
			int roughness = (roughnessData[row][col]>=0.5)? 1:0;
			float fel = felData[row][col];
			//cout << "row = " << row << " col= " << col << " sourceDir= " << (int)sourceDir << " floodProb= " << floodProb << endl;
			//tree.insert(row, col, (int)sourceDir, floodProb);
			int currentId = row * NCOLS + col;
			int parentId;
			switch (sourceDir) {
			case 1:
				parentId = currentId + 1;
				break;
			case 2:
				parentId = currentId + parameter.COLUMN + 1;
				break;
			case 4:
				parentId = currentId + parameter.COLUMN;
				break;
			case 8:
				parentId = currentId + parameter.COLUMN - 1;
				break;
			case 16:
				parentId = currentId - 1;
				break;
			case 32:
				parentId = currentId - parameter.COLUMN - 1;
				break;
			case 64:
				parentId = currentId - parameter.COLUMN;
				break;
			case 128:
				parentId = currentId - parameter.COLUMN + 1;
				break;
			default:
				continue;
			}
			int newCurrentId, newParentId;
			if (nodeIndexMap.find(currentId) == nodeIndexMap.end()) {
				nodeIndexMap.insert(make_pair(currentId, nodeIdCounter));
				newCurrentId = nodeIdCounter;
				data.allNodes.push_back(new Node(0, nodeIdCounter));
				data.allNodes[newCurrentId]->originalId = currentId;
				nodeIdCounter++;
			}
			else {
				newCurrentId = nodeIndexMap[currentId];
			}
			if (nodeIndexMap.find(parentId) == nodeIndexMap.end()) {
				nodeIndexMap.insert(make_pair(parentId, nodeIdCounter));
				//nodeIndexMap[parentId] = nodeIdCounter;
				newParentId = nodeIdCounter;
				data.allNodes.push_back(new Node(0, nodeIdCounter));
				data.allNodes[newParentId]->originalId = parentId;
				nodeIdCounter++;
			}
			else {
				newParentId = nodeIndexMap[parentId];
			}
			data.allNodes[newCurrentId]->parentsID.push_back(newParentId);   //source direction layer
			data.allNodes[newParentId]->childrenID.push_back(newCurrentId);  //source direction layer
			data.allNodes[newCurrentId]->bank = bank;                 //bank.tif
			data.allNodes[newCurrentId]->isObserved = (pits==1 || floodProb==0)? 0 : 1;           // no pits no na
			data.allNodes[newCurrentId]->isTree = tree;               // treecanopy.tif
			data.allNodes[newCurrentId]->isPits = pits;               // pits.tif
			data.allNodes[newCurrentId]->isNa = (floodProb==0)? 1:0;                 //na?  0-255 deltaresult.tif  0 is na
			data.allNodes[newCurrentId]->roughness = roughness;            // another input
			data.allNodes[newCurrentId]->cost = cost;                //cost.tif
			data.allNodes[newCurrentId]->p = floodProb/255.0;                   // delta
			data.allNodes[newCurrentId]->fel = fel;                 //field elevation layer
		}
	}
	//Added by Arpan for fist model
	//string line;
	//int nodeIdCounter = 0;
	//std::map<int, int> nodeIndexMap;
	//while (getline(sourceDirectionFile, line)) {
	//	//std::cout << line<<endl;
	//	vector<int> result;
	//	vector<float> fresult;
	//	stringstream s_stream(line);
	//	int counter = 0;
	//	while (s_stream.good()) {      //result 0-currentid, 1-sourceDirection 2-bank, 3- observation indicator 4-cost, 5-probability 
	//		string substr;
	//		getline(s_stream, substr, ',');
	//		if (counter < 8) {
	//			result.push_back(stoi(substr));
	//			counter++;
	//		}
	//		else {
	//			fresult.push_back(stof(substr));
	//			counter++;
	//		}
	//	}
	//	int currentId = result[0];
	//	int parentId;
	//	switch (result[1]) {
	//	case 1:
	//		parentId = currentId + 1;
	//		break;
	//	case 2:
	//		parentId = currentId + parameter.COLUMN + 1;
	//		break;
	//	case 4:
	//		parentId = currentId + parameter.COLUMN;
	//		break;
	//	case 8:
	//		parentId = currentId + parameter.COLUMN - 1;
	//		break;
	//	case 16:
	//		parentId = currentId - 1;
	//		break;
	//	case 32:
	//		parentId = currentId - parameter.COLUMN - 1;
	//		break;
	//	case 64:
	//		parentId = currentId - parameter.COLUMN;
	//		break;
	//	case 128:
	//		parentId = currentId - parameter.COLUMN + 1;
	//		break;
	//	default:
	//		std::cout << "Error in source direction" << std::endl;

	//	}
	//	int newCurrentId, newParentId;
	//	if (nodeIndexMap.find(currentId) == nodeIndexMap.end()) {
	//		nodeIndexMap.insert(make_pair(currentId, nodeIdCounter));
	//		//nodeIndexMap[currentId] = nodeIdCounter;
	//		newCurrentId = nodeIdCounter;
	//		data.allNodes.push_back(new Node(0, nodeIdCounter));
	//		data.allNodes[newCurrentId]->originalId = currentId;
	//		nodeIdCounter++;
	//	}
	//	else {
	//		newCurrentId = nodeIndexMap[currentId];
	//	}
	//	if (nodeIndexMap.find(parentId) == nodeIndexMap.end()) {
	//		nodeIndexMap.insert(make_pair(parentId, nodeIdCounter));
	//		//nodeIndexMap[parentId] = nodeIdCounter;
	//		newParentId = nodeIdCounter;
	//		data.allNodes.push_back(new Node(0, nodeIdCounter));
	//		data.allNodes[newParentId]->originalId = parentId;
	//		nodeIdCounter++;
	//	}
	//	else {
	//		newParentId = nodeIndexMap[parentId];
	//	}
	//	data.allNodes[newCurrentId]->parentsID.push_back(newParentId);   //source direction layer
	//	data.allNodes[newParentId]->childrenID.push_back(newCurrentId);  //source direction layer
	//	data.allNodes[newCurrentId]->bank = result[2];                 //bank.tif
	//	data.allNodes[newCurrentId]->isObserved = result[3];           // no pits no trees no na
	//	data.allNodes[newCurrentId]->isTree = result[4];               // treecanopy.tif
	//	data.allNodes[newCurrentId]->isPits = result[5];               // pits.tif
	//	data.allNodes[newCurrentId]->isNa = result[6];                 //na?  0 -255 deltaresult.tif  0 is na
	//	data.allNodes[newCurrentId]->roughness = result[7];            // another input
	//	data.allNodes[newCurrentId]->cost = fresult[0];                //cost.tif
	//	data.allNodes[newCurrentId]->p = fresult[1];                   // delta
	//	data.allNodes[newCurrentId]->fel = fresult[2];                 //field elevation layer
	//}
	//adjust reach ids with new index structure 
	//Note to Saramsha: automatically handled by your case 
	for (int i = 0; i < data.reach_ids.size(); i++) {
		int currentreachId = data.reach_ids[i];
		int row = (int)(currentreachId / parameter.COLUMN);
		int col = currentreachId % parameter.COLUMN;
		float fel = felData[row][col];
		data.reach_fel.push_back(fel);
		int newCurrentId = -1;
		if (nodeIndexMap.find(currentreachId) == nodeIndexMap.end()) {
			nodeIndexMap.insert(make_pair(currentreachId, nodeIdCounter));
			//nodeIndexMap[currentId] = nodeIdCounter;
			newCurrentId = nodeIdCounter;
			data.allNodes.push_back(new Node(0, nodeIdCounter));
			data.allNodes[newCurrentId]->originalId = currentreachId;
			data.AdjustedReachNodes.push_back(newCurrentId);
			data.allNodes[newCurrentId]->bank = 0;
			data.allNodes[newCurrentId]->cost = 0.0;
			data.allNodes[newCurrentId]->p = 1.0;
			data.allNodes[newCurrentId]->isObserved = 0; // to indicate reach nodes
			data.allNodes[newCurrentId]->fel = fel;
			nodeIdCounter++;

		}
		else {
			newCurrentId = nodeIndexMap[currentreachId];
			data.AdjustedReachNodes.push_back(newCurrentId);
			data.allNodes[newCurrentId]->bank = 0;
			data.allNodes[newCurrentId]->cost = 0.0;
			data.allNodes[newCurrentId]->p = 1.0;
			data.allNodes[newCurrentId]->isObserved = 0; // to indicate reach nodes
			data.allNodes[newCurrentId]->fel = fel;
		}
	}

	data.inferedmaxCostLeft.resize(data.AdjustedReachNodes.size(), -1);
	data.inferedmaxCostRight.resize(data.AdjustedReachNodes.size(), -1);

	//get root node for each tree
	data.leftbfsRootNodes.resize(data.AdjustedReachNodes.size(), -1);  // leaf nodes with no children
	data.rightbfsRootNodes.resize(data.AdjustedReachNodes.size(), -1); //leaf nodes with no children
	for (int i = 0; i < data.AdjustedReachNodes.size(); i++) {
		int nodeId = data.AdjustedReachNodes[i];
		//finding root in left tree
		int leftnode = -1;
		int rightnode = -1;
		for (int j = 0; j < data.allNodes[nodeId]->childrenID.size(); j++) {
			int cid = data.allNodes[nodeId]->childrenID[j];
			if (data.allNodes[cid]->bank == 1) {
				leftnode = cid;
				break;
			}

		}
		for (int j = 0; j < data.allNodes[nodeId]->childrenID.size(); j++) {
			int cid = data.allNodes[nodeId]->childrenID[j];
			if (data.allNodes[cid]->bank == 2) {
				rightnode = cid;
				break;
			}
		}
		if (leftnode != -1) {
			while (data.allNodes[leftnode]->childrenID.size() != 0) {
				leftnode = data.allNodes[leftnode]->childrenID[0];
			}
		}
		data.leftbfsRootNodes[i] = leftnode;
		if (rightnode != -1) {
			while (data.allNodes[rightnode]->childrenID.size() != 0) {
				rightnode = data.allNodes[rightnode]->childrenID[0];
			}
		}
		data.rightbfsRootNodes[i] = rightnode;

	}

	//get bfs order for each tree

	vector<int> bfsVisited;
	bfsVisited.resize(data.allNodes.size(), 0);
	for (int i = 0; i < data.AdjustedReachNodes.size(); i++) {
		if (data.leftbfsRootNodes[i] == -1) {
			data.leftbfsOrder.push_back({});
		}
		else {
			data.leftbfsOrder.push_back(getBFSOrder(data.leftbfsRootNodes[i], bfsVisited, 1));
		}
	}
	std::fill(bfsVisited.begin(), bfsVisited.end(), 0);
	for (int i = 0; i < data.AdjustedReachNodes.size(); i++) {
		if (data.rightbfsRootNodes[i] == -1) {
			data.rightbfsOrder.push_back({});
		}
		else {
			data.rightbfsOrder.push_back(getBFSOrder(data.rightbfsRootNodes[i], bfsVisited, 2));
		}
	}

	data.hasObservedPixelsLeft.resize(data.leftbfsOrder.size(), 0);
	data.hasObservedPixelsRight.resize(data.rightbfsOrder.size(), 0);
	//std::cout << endl << "Left: ";
	for (int i = 0; i < data.leftbfsOrder.size(); i++) {
		for (int treeIndex = 0; treeIndex < data.leftbfsOrder[i].size(); treeIndex++) {
			int pixelId = data.leftbfsOrder[i][treeIndex];
			if (data.allNodes[pixelId]->isObserved == 1) {
				data.hasObservedPixelsLeft[i] = 1;
				//std::cout << i << " ";
				break;
			}
		}
	}
	//std::cout << endl << "Right: ";
	for (int i = 0; i < data.rightbfsOrder.size(); i++) {
		for (int treeIndex = 0; treeIndex < data.rightbfsOrder[i].size(); treeIndex++) {
			int pixelId = data.rightbfsOrder[i][treeIndex];
			if (data.allNodes[pixelId]->isObserved == 1) {
				data.hasObservedPixelsRight[i] = 1;
				//std::cout << i << " ";
				break;
			}
		}
	}

	for (int i = 0; i < data.AdjustedReachNodes.size(); i++) {
		int nodeId = data.AdjustedReachNodes[i];
		if (data.allNodes[nodeId]->originalId != data.reach_ids[i]) {
			cout <<"Error " <<data.allNodes[nodeId]->originalId << "!= " << data.reach_ids[i];
			exit(0);
		}
	}
	//ifstream ReachFELFile(CTInputLocation + CTReachFEl);
	//string fline;
	//int index=0;
	//while (getline(ReachFELFile, fline)) {
	//	int rNode = data.AdjustedReachNodes[index];
	//	std::cout.precision(16);
	//	double fel = stod(fline);
	//	data.allNodes[rNode]->fel=fel;
	//	index++;
	//}
	//ReachFELFile.close();
	data.highestCostLeft.resize(data.AdjustedReachNodes.size(), MAXCOST);
	data.highestCostRight.resize(data.AdjustedReachNodes.size(), MAXCOST);
	for (int i = 0; i < data.leftbfsOrder.size(); i++) {
		for (int treeIndex = 0; treeIndex < data.leftbfsOrder[i].size(); treeIndex++) {
			int pixelId = data.leftbfsOrder[i][treeIndex];
			if (data.allNodes[pixelId]->cost > data.highestCostLeft[i]) {
				data.highestCostLeft[i] = data.allNodes[pixelId]->cost;
			}
		}
	}
	for (int i = 0; i < data.rightbfsOrder.size(); i++) {
		for (int treeIndex = 0; treeIndex < data.rightbfsOrder[i].size(); treeIndex++) {
			int pixelId = data.rightbfsOrder[i][treeIndex];
			if (data.allNodes[pixelId]->cost > data.highestCostRight[i]) {
				data.highestCostRight[i] = data.allNodes[pixelId]->cost;
			}
		}
	}
	parameter.allPixelSize = data.allNodes.size();
	parameter.elnPzn_xn.resize(parameter.allPixelSize * cNum, eln(0.5));

	for (int i = 0; i < parameter.allPixelSize; i++) {
		if (data.allNodes[i]->isObserved == -1) {
			parameter.elnPzn_xn[i * cNum] = eln(0.5);
			parameter.elnPzn_xn[i * cNum + 1] = eln(0.5);
		}
		else if (data.allNodes[i]->isObserved == 0) {  //stream nodes
			data.allNodes[i]->p = 0.999;
			parameter.elnPzn_xn[i * cNum] = eln(0.001);
			parameter.elnPzn_xn[i * cNum + 1] = eln(0.999);
		}
		else {
			float probability = data.allNodes[i]->p;
			if (probability < 0 || probability>1) {
				std::cout << "i= " << i << " prob erorr :" << probability << endl;
			}
			parameter.elnPzn_xn[i * cNum] = eln(1 - probability);
			parameter.elnPzn_xn[i * cNum + 1] = eln(probability);
		}
	}

	for (int i = 0; i < data.leftbfsRootNodes.size(); i++) {
		if (data.leftbfsRootNodes[i] != -1) {
			data.leftbfsRootNodes_orgId.push_back(data.allNodes[data.leftbfsRootNodes[i]]->originalId);
		}
		else {
			data.leftbfsRootNodes_orgId.push_back(-1);
		}
		if (data.rightbfsRootNodes[i] != -1) {
			data.rightbfsRootNodes_orgId.push_back(data.allNodes[data.rightbfsRootNodes[i]]->originalId);
		}
		else {
			data.rightbfsRootNodes_orgId.push_back(-1);
		}

	}
	//// Test: remove later
	//getIds();
	//getOrgIds();
	//getOriginIdBanks();
	//getOriginIdLeftBanks();
	//getOriginIdRightBanks();
	//getRegionNodeCount();
	//getLeftRegionNodeCount();
	//getRightRegionNodeCount();
	////end Test

	//convert parameter Pi, M, Epsilon to log form
	parameter.Pi = eln(parameter.Pi);
	parameter.Epsilon = eln(parameter.Epsilon); //check if already eln form?

	this->UpdateTransProb();
	string reachId = CTPara.substr(0, 2);
	parameter.reachId = reachId;
	string EffectiveBrach = "CL";
	if (EB == 1) {
		EffectiveBrach = "MC";
	}
	else if (EB == 2) {
		EffectiveBrach = "BC";
	}
	parameter.fname = reachId + "_pr_DV3_Top" + to_string((int)(parameter.cutoff_percentage * 100)) + "_" + EffectiveBrach + "_BN" + to_string(BOUNDARY_NODES_OBSERVED);


	std::cout << "Inference Started.." << endl;
	auto start = std::chrono::system_clock::now();
	inference();
	//sanityChecker();
	if (EFFECTIVE_BRANCH_VIZ == 1) {
		getOriginIdBanks_effectiveBranches();
	}
	for (int i = 0; i < data.AdjustedReachNodes.size(); i++) {
		data.inferedcostLeft_afterInference.push_back(data.inferedmaxCostLeft[i]);
		data.inferedcostRight_afterInference.push_back(data.inferedmaxCostRight[i]);
	}

	interpolate();

	prediction();
	auto end = std::chrono::system_clock::now();
	auto elapsed_seconds = end - start;

	std::cout << "Inference Finished. Duration: " << elapsed_seconds.count() << endl << endl;
	output();

}
void cFlood::prediction() {
	//mappredictions.resize(parameter.orgPixelSize, -1);
	std::fill(mappredictions.begin(), mappredictions.end(), -1);
	//for left trees 
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		//if (data.hasObservedPixelsLeft[leftOrder] && data.leftbfsOrder[leftOrder].size()>=PIXELLIMT) {
		//	continue;
		//}
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				if (data.inferedmaxCostLeft[leftOrder] == -1 && data.leftbfsOrder[leftOrder].size() < PIXELLIMT) {  //not interpolated
					if (data.allNodes[nodid]->isNa == 0)
						mappredictions[data.allNodes[nodid]->originalId] = 1;
				}
				else {
					if (data.allNodes[nodid]->cost <= data.inferedmaxCostLeft[leftOrder] && data.allNodes[nodid]->isNa==0) {
						mappredictions[data.allNodes[nodid]->originalId] = 1;
					}
					else if(data.allNodes[nodid]->isNa ==0) {
						mappredictions[data.allNodes[nodid]->originalId] = 0;
					}
				}
			}
		}
	}

	//for right trees 
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		//if (data.hasObservedPixelsRight[rightOrder] && data.rightbfsOrder[rightOrder].size()>= PIXELLIMT) {
		//	continue;
		//}
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				if (data.inferedmaxCostRight[rightOrder] == -1 && data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {  //not interpolated
					if(data.allNodes[nodid]->isNa == 0)
						mappredictions[data.allNodes[nodid]->originalId] = 1;
				}
				else {
					if (data.allNodes[nodid]->cost <= data.inferedmaxCostRight[rightOrder] && data.allNodes[nodid]->isNa == 0) {
						mappredictions[data.allNodes[nodid]->originalId] = 1;
					}
					else if(data.allNodes[nodid]->isNa == 0){
						mappredictions[data.allNodes[nodid]->originalId] = 0;
					}
				}
			}
		}
	}
	for (int i = 0; i < data.reach_ids.size(); i++) {
		mappredictions[data.reach_ids[i]] = 1;  //reach ids are the lowest point in the river
	}
}

void cFlood::interpolate() {
	//profile table before interpolation
	ofstream profiletable;
	data.combinedCost.resize(data.reach_ids.size(), 0);
	data.avgCost.resize(data.reach_ids.size(), 0);
	//find the regions with infered cost values (i.e values that are not -1)
	vector<int> stops;
	for (int i = 0; i < data.reach_ids.size(); i++) {
		if (data.inferedmaxCostLeft[i] > 0 || data.inferedmaxCostRight[i] > 0) {
			stops.push_back(i);
		}
	}
	//calculating combined cost
	//case 1: from first index to first stop
	float value = (data.inferedmaxCostLeft[stops[0]] > 0) ? value = data.inferedmaxCostLeft[stops[0]] : data.inferedmaxCostRight[stops[0]];
	for (int i = 0; i <= stops[0]; i++) {
		data.combinedCost[i] = value;
	}
	//case 2: from last infered cost to last reach ids
	int lastIndex = (stops.size() - 1);
	value = (data.inferedmaxCostLeft[stops[lastIndex]] > 0) ? data.inferedmaxCostLeft[stops[lastIndex]] : data.inferedmaxCostRight[stops[lastIndex]];
	for (int i = stops[lastIndex]; i < data.reach_ids.size(); i++) {
		data.combinedCost[i] = value;
	}
	//case 3: intermediate
	for (int i = 0; i < stops.size() - 1; i++) {
		//cout<<"i = "<<i<<endl;
		int first = stops[i];
		int second = stops[(i + 1)];
		int diff = second - first;
		//cout << "i = " << i << " first= " <<first<<" second= "<<second <<" diff= "<<diff<< endl;
		//cout << " data.inferedmaxCostLeft[stops[second]] = " << data.inferedmaxCostLeft[stops[second]] << " data.inferedmaxCostRight[stops[second]]= " << data.inferedmaxCostRight[stops[second]]<< endl;
		float firstValue = (data.inferedmaxCostLeft[first] > 0) ? data.inferedmaxCostLeft[first] : data.inferedmaxCostRight[first];
		float secondValue = (data.inferedmaxCostLeft[second] > 0) ? data.inferedmaxCostLeft[second] : data.inferedmaxCostRight[second];
		//cout << "i = " << i << " first= " << first << " second= " << second << " diff= " << diff << " firstvalue= " << firstValue << " secondValue= " << secondValue << endl;
		float change = (secondValue - firstValue) / diff;
		data.combinedCost[first] = firstValue;
		data.combinedCost[second] = secondValue;
		for (int j = first + 1; j < second; j++) {
			data.combinedCost[j] = data.combinedCost[(j - 1)] + change;
		}
	}

	profiletable.open(CTOutputLocation + "ProfileTables\\" + parameter.reachId + "_ProfileTable_preInterpolation.csv");
	//profiletable.open(CTOutputLocation + "ProfileTables\\" + parameter.fname + "_ProfileTable_preInterpolation.csv");
	profiletable << "SourceId" << "," << "fel" << "," << "Distance Covered" << "," << "Flow Accumulation" << "," << "Left Node Count" << "," << "Left Max Cost" << "," << "Left Infered Cost" << ","
		<< "Observation Indicator" << "," << "Right Node Count" << "," << "Right Max Cost" << "," << "Right Infered Cost" << ","
		<< "Observation Indicator" << endl;
	for (int index = 0; index < data.AdjustedReachNodes.size(); index++) {
		profiletable << data.reach_ids[index] << ","
			<< data.reach_fel[index] << ","
			<< data.distance_covered[index] << ","
			<< data.flow_accumulation[index] << ","
			//<< data.allNodes[data.AdjustedReachNodes[index]]->fel<<","
			<< data.leftbfsOrder[index].size() << ","
			<< data.highestCostLeft[index] << ","
			<< data.inferedmaxCostLeft[index] << ","
			<< data.hasObservedPixelsLeft[index] << ","
			<< data.rightbfsOrder[index].size() << ","
			<< data.highestCostRight[index] << ","
			<< data.inferedmaxCostRight[index] << ","
			<< data.hasObservedPixelsRight[index] << endl;
	}
	profiletable.close();
	int current = 0;
	//for left;
	while (current < data.AdjustedReachNodes.size()) {
		if (data.inferedmaxCostLeft[current] == -1 && current == 0) {
			//find the first reach node with non -1 max cost value
			int index = -1;
			for (int j = 1; j < data.AdjustedReachNodes.size(); j++) {
				if (data.inferedmaxCostLeft[j] != -1) {
					index = j;
					break;
				}
			}
			if (index == -1) {
				break;
			}
			double value = data.inferedmaxCostLeft[index];
			for (int i = 0; i < index; i++) {
				data.inferedmaxCostLeft[i] = value;
			}
			current = index;
		}
		else if (data.inferedmaxCostLeft[current] != -1) {
			//two cases
				//case 1: there are n points in between next reach that has cost value 
				//case 2: there is no next point
			//find index of next reach node that has cost value
			int index = -1;
			int count = 0;
			double value = data.inferedmaxCostLeft[current];
			for (int j = current + 1; j < data.AdjustedReachNodes.size(); j++) {
				if (data.inferedmaxCostLeft[j] != -1) {
					index = j;
					break;
				}
				count++;
			}
			if (index == -1) {// case 2
				for (int i = current + 1; i < data.AdjustedReachNodes.size(); i++) {
					data.inferedmaxCostLeft[i] = value;
				}
				current = data.AdjustedReachNodes.size();
				break;
			}
			else if (count == 0 && index == current + 1) {
				current = index;
			}
			else {
				double interval = (data.inferedmaxCostLeft[index] - value) / count;
				for (int i = current + 1; i < index; i++) {
					data.inferedmaxCostLeft[i] = data.inferedmaxCostLeft[(i - 1)] + interval;
				}
				current = index;
			}
		}

	}

	//for right bank.
	current = 0;
	while (current < data.AdjustedReachNodes.size()) {
		if (data.inferedmaxCostRight[current] == -1 && current == 0) {
			//find the first reach node with non -1 max cost value
			int index = -1;
			for (int j = 1; j < data.AdjustedReachNodes.size(); j++) {
				if (data.inferedmaxCostRight[j] != -1) {
					index = j;
					break;
				}
			}
			if (index == -1) {
				break;
			}
			double value = data.inferedmaxCostRight[index];
			for (int i = 0; i < index; i++) {
				data.inferedmaxCostRight[i] = value;
			}
			current = index;
		}
		else if (data.inferedmaxCostRight[current] != -1) {
			//two cases
				//case 1: there are n points in between next reach that has cost value 
				//case 2: there is no next point
			//find index of next reach node that has cost value
			int index = -1;
			int count = 0;
			double value = data.inferedmaxCostRight[current];
			for (int j = current + 1; j < data.AdjustedReachNodes.size(); j++) {
				if (data.inferedmaxCostRight[j] != -1) {
					index = j;
					break;
				}
				count++;
			}
			if (index == -1) {// case 2
				for (int i = current + 1; i < data.AdjustedReachNodes.size(); i++) {
					data.inferedmaxCostRight[i] = value;
				}
				current = data.AdjustedReachNodes.size();
				break;
			}
			else if (count == 0 && index == current + 1) {
				current = index;
			}
			else {
				double interval = (data.inferedmaxCostRight[index] - value) / count;
				for (int i = current + 1; i < index; i++) {
					data.inferedmaxCostRight[i] = data.inferedmaxCostRight[(i - 1)] + interval;
				}
				current = index;
			}
		}

	}

}

void cFlood::removeLink(vector<int>& v, int removeID) {
	v.erase(std::find(v.begin(), v.end(), removeID));
}

void cFlood::inference() {
	vector<int> inferVisited(parameter.allPixelSize, 0);
	for (int i = 0; i < parameter.allPixelSize; i++) {
		data.allNodes[i]->correspondingNeighbour.clear();
		data.allNodes[i]->correspondingNeighbourClassOne.clear();
		data.allNodes[i]->correspondingNeighbourClassZero.clear();
	}
	//data.allNodes.correspondingNeighboursClass.resize(parameter.allPixelSize * cNum);
	//infer.correspondingNeighbours.resize(parameter.allPixelSize);
	//for left
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (!data.hasObservedPixelsLeft[leftOrder] || data.leftbfsOrder[leftOrder].size() < PIXELLIMT) {
			continue;
		}
		int bfsTraversalOrderSize = (int)data.leftbfsOrder[leftOrder].size();
		//int bfsTraversalOrderSize = (int)data.bfsTraversalOrder.size();
		for (int node = bfsTraversalOrderSize - 1; node >= 0; node--) {
			int cur_node_id = data.leftbfsOrder[leftOrder][node];
			//std::cout << "cur id = " << cur_node_id << endl;
			//what to do in case of reach node?
			//if (cur_node_id == data.AdjustedReachNodes[leftOrder]) {
			//
			//}
			//int cur_node_id = data.bfsTraversalOrder[node];
			vector<int> leftbankchildrenID;
			for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
				int child = data.allNodes[cur_node_id]->childrenID[c];
				if (data.allNodes[child]->bank == 1) {
					leftbankchildrenID.push_back(child);
				}
			}
			//data.allNodes[cur_node_id]->fi_ChildList.resize(data.allNodes[cur_node_id]->childrenID.size() * cNum, 0);
			data.allNodes[cur_node_id]->fi_ChildList.resize(leftbankchildrenID.size() * cNum, 0);
			for (int cls = 0; cls < cNum; cls++) {
				data.allNodes[cur_node_id]->fi_parent[cls] = 0;
				data.allNodes[cur_node_id]->fo[cls] = 0;
			}

			//first figure out which neighbor fmessage passes to from current node pass n->? foNode;
			//idea: In bfs traversal order leave to root, check if next the node in bfs order is parent or child of the current node (should be child or parent of the current node)
			int foNode = -1;
			bool foNode_isChild = false;
			for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {  //check in parent list if found respective parent node is foNode
				int pid = data.allNodes[cur_node_id]->parentsID[p];
				if (!inferVisited[pid]) {
					foNode = pid;
					break;
				}
			}
			if (foNode == -1) {
				//for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
				for (int c = 0; c < leftbankchildrenID.size(); c++) {
					int cid = data.allNodes[cur_node_id]->childrenID[c];
					if (!inferVisited[cid]) {
						foNode = cid;
						foNode_isChild = true;
						break;
					}
				}
			}
			data.allNodes[cur_node_id]->foNode = foNode;
			data.allNodes[cur_node_id]->foNode_ischild = foNode_isChild;
			if (cur_node_id == data.leftbfsRootNodes[leftOrder] && leftbankchildrenID.size() == 0) { //need to verify this changed && to || for IRONFIST project
				foNode_isChild = true;
			}

			//incoming message from visited child
			if (leftbankchildrenID.size() > 0) {

				for (int c = 0; c < leftbankchildrenID.size(); c++) {
					int child_id = leftbankchildrenID[c];

					if (child_id == foNode) {
						continue;
					}
					data.allNodes[cur_node_id]->correspondingNeighbour.push_back(child_id);
					for (int p = 0; p < data.allNodes[child_id]->parentsID.size(); p++) {
						int pid = data.allNodes[child_id]->parentsID[p];
						if (pid != cur_node_id) {
							data.allNodes[cur_node_id]->correspondingNeighbour.push_back(pid);
						}

					}
					vector<int> parentOfChildExcept_currentNode;
					for (int en = 0; en < data.allNodes[child_id]->parentsID.size(); en++) {
						if (data.allNodes[child_id]->parentsID[en] != cur_node_id) {
							parentOfChildExcept_currentNode.push_back(data.allNodes[child_id]->parentsID[en]);
						}

					}
					for (int cls = 0; cls < cNum; cls++) {  //cls represents current node class
															//double sumAccumulator = eln(0);   //should be 0 since we are summing it up//eln(1);//need to confirm
						double max = eln(0);
						vector<int> maxCorrespondingNeighbour;
						for (int c_cls = 0; c_cls < cNum; c_cls++) { //c_cls reperesnets child class label   Yc
							int max_bitCount = 1 << parentOfChildExcept_currentNode.size();
							for (int bitCount = 0; bitCount < max_bitCount; bitCount++) { //summation for each parent and child class label(given by c_cls)
								double productAccumulator = data.allNodes[child_id]->fo[c_cls];  //product with fo(c)
								vector<int>neighbourClass;
								neighbourClass.push_back(c_cls);
								int parentClsProd = 1; //p(c), product of parent classes for child c
								for (int p = 0; p < parentOfChildExcept_currentNode.size(); p++) {//calculating Product(fo(p)) for all parent of current child except the current node
									int pid = parentOfChildExcept_currentNode[p];
									int parentClsValue = (bitCount >> p) & 1;
									parentClsProd *= parentClsValue;
									neighbourClass.push_back(parentClsValue);
									productAccumulator = elnproduct(productAccumulator, data.allNodes[pid]->fo[parentClsValue]);  //calculating Product(fo(p)) for all parent of current child except the current node
								}
								//multiplying P(Yc|Ypc)
								parentClsProd *= cls;
								productAccumulator = elnproduct(productAccumulator, parameter.elnPz_zpn[c_cls][parentClsProd]);
								if (max < productAccumulator) {
									max = productAccumulator;
									maxCorrespondingNeighbour = neighbourClass;
								}
							}
						}
						data.allNodes[cur_node_id]->fi_ChildList[(c * cNum + cls)] = max;
						if (cls == 0) {
							for (int t = 0; t < maxCorrespondingNeighbour.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassZero.push_back(maxCorrespondingNeighbour[t]);
							}
						}
						else {
							for (int t = 0; t < maxCorrespondingNeighbour.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassOne.push_back(maxCorrespondingNeighbour[t]);
							}
						}
					}
				}
			}

			if (foNode_isChild) {  //means the current node has all visited parents
				if (data.allNodes[cur_node_id]->parentsID.size() == 0) {
					for (int cls = 0; cls < cNum; cls++) {
						data.allNodes[cur_node_id]->fi_parent[cls] = parameter.elnPz[cls];
					}
				}
				else {
					for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {
						int pid = data.allNodes[cur_node_id]->parentsID[p];
						data.allNodes[cur_node_id]->correspondingNeighbour.push_back(pid);
					}
					for (int cls = 0; cls < cNum; cls++) {
						double max = eln(0);
						vector<int> maxNeighbourClass;
						int max_bitCount = 1 << data.allNodes[cur_node_id]->parentsID.size();
						for (int bitCount = 0; bitCount < max_bitCount; bitCount++) { //summation for each parent class label
							vector<int> parentClass;
							double productAccumulator = eln(1);
							int parentClsProd = 1;
							for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {
								int pid = data.allNodes[cur_node_id]->parentsID[p];
								int parentClsValue = (bitCount >> p) & 1;
								parentClass.push_back(parentClsValue);
								parentClsProd *= parentClsValue;
								productAccumulator = elnproduct(productAccumulator, data.allNodes[pid]->fo[parentClsValue]);  //calculating Product(fo(p)) for all parent of current child except the current node
							}
							productAccumulator = elnproduct(productAccumulator, parameter.elnPz_zpn[cls][parentClsProd]);
							if (max < productAccumulator) {
								max = productAccumulator;
								maxNeighbourClass = parentClass;
							}
							//sumAccumulator = elnsum(sumAccumulator, productAccumulator);
						}
						data.allNodes[cur_node_id]->fi_parent[cls] = max;
						if (cls == 0) {
							for (int t = 0; t < maxNeighbourClass.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassZero.push_back(maxNeighbourClass[t]);
							}
						}
						else {
							for (int t = 0; t < maxNeighbourClass.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassOne.push_back(maxNeighbourClass[t]);
							}
						}
					}
				}

				//calulating fo
				for (int cls = 0; cls < cNum; cls++) { //cls represents class of the current node
					double productAccumulator = eln(1);
					for (int c = 0; c < leftbankchildrenID.size(); c++) {
						int child_id = leftbankchildrenID[c];
						if (child_id == foNode) {
							continue;
						}
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_ChildList[c * cNum + cls]); //multiplying with al the child fi except the outgoing child
					}
					productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_parent[cls]);  // multiplying with fi(n)_parent
					productAccumulator = elnproduct(productAccumulator, parameter.elnPzn_xn[(cur_node_id * cNum + cls)]);
					productAccumulator = elnproduct(productAccumulator, -1 * parameter.elnPz[cls]);
					data.allNodes[cur_node_id]->fo[cls] = productAccumulator;
				}

			}

			else {  //message pass n-> parent there is no fi(n)_parent   //computes for root node as well
					//calulating fo
				for (int cls = 0; cls < cNum; cls++) { //cls represents class of the current node
					double productAccumulator = eln(1);
					for (int c = 0; c < leftbankchildrenID.size(); c++) {
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_ChildList[(c * cNum + cls)]); //multiplying with al the child fi except the outgoing child
					}
					//productAccumulator = elnproduct(productAccumulator, parameter.elnPxn_zn[cur_node_id*cNum + cls]);
					productAccumulator = elnproduct(productAccumulator, parameter.elnPzn_xn[(cur_node_id * cNum + cls)]);
					productAccumulator = elnproduct(productAccumulator, -1 * parameter.elnPz[cls]);
					data.allNodes[cur_node_id]->fo[cls] = productAccumulator;
				}
			}

			inferVisited[cur_node_id] = 1;
		}
	}
	//need to run map prediction for left before runing inference
	updateMapPrediction_left();
	//
	//	//for right
	inferVisited.clear();
	inferVisited.resize(parameter.allPixelSize, 0);
	for (int i = 0; i < parameter.allPixelSize; i++) {
		data.allNodes[i]->correspondingNeighbour.clear();
		data.allNodes[i]->correspondingNeighbourClassOne.clear();
		data.allNodes[i]->correspondingNeighbourClassZero.clear();
	}
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (!data.hasObservedPixelsRight[rightOrder] || data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
			continue;
		}
		int bfsTraversalOrderSize = (int)data.rightbfsOrder[rightOrder].size();
		//int bfsTraversalOrderSize = (int)data.bfsTraversalOrder.size();
		for (int node = bfsTraversalOrderSize - 1; node >= 0; node--) {
			int cur_node_id = data.rightbfsOrder[rightOrder][node];
			//int cur_node_id = data.bfsTraversalOrder[node];
			vector<int> rightbankchildrenID;
			for (int c = 0; c < data.allNodes[cur_node_id]->childrenID.size(); c++) {
				int child = data.allNodes[cur_node_id]->childrenID[c];
				if (data.allNodes[child]->bank == 2) {
					rightbankchildrenID.push_back(child);
				}
			}
			//			data.allNodes[cur_node_id]->fi_ChildList.resize(data.allNodes[cur_node_id]->childrenID.size()* cNum, 0);
			data.allNodes[cur_node_id]->fi_ChildList.resize(rightbankchildrenID.size() * cNum, 0);
			for (int cls = 0; cls < cNum; cls++) {
				data.allNodes[cur_node_id]->fi_parent[cls] = 0;
				data.allNodes[cur_node_id]->fo[cls] = 0;
			}

			//first figure out which neighbor fmessage passes to from current node pass n->? foNode;
			//idea: In bfs traversal order leave to root, check if next the node in bfs order is parent or child of the current node (should be child or parent of the current node)
			int foNode = -1;
			bool foNode_isChild = false;
			for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {  //check in parent list if found respective parent node is foNode
				int pid = data.allNodes[cur_node_id]->parentsID[p];
				if (!inferVisited[pid]) {
					foNode = pid;
					break;
				}
			}
			if (foNode == -1) {
				for (int c = 0; c < rightbankchildrenID.size(); c++) {
					int cid = rightbankchildrenID[c];
					if (!inferVisited[cid]) {
						foNode = cid;
						foNode_isChild = true;
						break;
					}
				}
			}
			data.allNodes[cur_node_id]->foNode = foNode;
			data.allNodes[cur_node_id]->foNode_ischild = foNode_isChild;
			if (cur_node_id == data.rightbfsRootNodes[rightOrder] && rightbankchildrenID.size() == 0) { //need to verify this changed && to || for IRONFIST project
				foNode_isChild = true;
			}

			//incoming message from visited child
			if (data.allNodes[cur_node_id]->childrenID.size() > 0) {

				for (int c = 0; c < rightbankchildrenID.size(); c++) {
					int child_id = rightbankchildrenID[c];

					if (child_id == foNode) {
						continue;
					}
					data.allNodes[cur_node_id]->correspondingNeighbour.push_back(child_id);
					for (int p = 0; p < data.allNodes[child_id]->parentsID.size(); p++) {
						int pid = data.allNodes[child_id]->parentsID[p];
						if (pid != cur_node_id) {
							data.allNodes[cur_node_id]->correspondingNeighbour.push_back(pid);
						}

					}
					vector<int> parentOfChildExcept_currentNode;
					for (int en = 0; en < data.allNodes[child_id]->parentsID.size(); en++) {
						if (data.allNodes[child_id]->parentsID[en] != cur_node_id) {
							parentOfChildExcept_currentNode.push_back(data.allNodes[child_id]->parentsID[en]);
						}

					}
					for (int cls = 0; cls < cNum; cls++) {  //cls represents current node class
															//double sumAccumulator = eln(0);   //should be 0 since we are summing it up//eln(1);//need to confirm
						double max = eln(0);
						vector<int> maxCorrespondingNeighbour;
						for (int c_cls = 0; c_cls < cNum; c_cls++) { //c_cls reperesnets child class label   Yc
							int max_bitCount = 1 << parentOfChildExcept_currentNode.size();
							for (int bitCount = 0; bitCount < max_bitCount; bitCount++) { //summation for each parent and child class label(given by c_cls)
								double productAccumulator = data.allNodes[child_id]->fo[c_cls];  //product with fo(c)
								vector<int>neighbourClass;
								neighbourClass.push_back(c_cls);
								int parentClsProd = 1; //p(c), product of parent classes for child c
								for (int p = 0; p < parentOfChildExcept_currentNode.size(); p++) {//calculating Product(fo(p)) for all parent of current child except the current node
									int pid = parentOfChildExcept_currentNode[p];
									int parentClsValue = (bitCount >> p) & 1;
									parentClsProd *= parentClsValue;
									neighbourClass.push_back(parentClsValue);
									productAccumulator = elnproduct(productAccumulator, data.allNodes[pid]->fo[parentClsValue]);  //calculating Product(fo(p)) for all parent of current child except the current node
								}
								//multiplying P(Yc|Ypc)
								parentClsProd *= cls;
								productAccumulator = elnproduct(productAccumulator, parameter.elnPz_zpn[c_cls][parentClsProd]);
								if (max < productAccumulator) {
									max = productAccumulator;
									maxCorrespondingNeighbour = neighbourClass;
								}
							}
						}
						data.allNodes[cur_node_id]->fi_ChildList[(c * cNum + cls)] = max;
						if (cls == 0) {
							for (int t = 0; t < maxCorrespondingNeighbour.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassZero.push_back(maxCorrespondingNeighbour[t]);
							}
						}
						else {
							for (int t = 0; t < maxCorrespondingNeighbour.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassOne.push_back(maxCorrespondingNeighbour[t]);
							}
						}
					}
				}
			}

			if (foNode_isChild) {  //means the current node has all visited parents
				if (data.allNodes[cur_node_id]->parentsID.size() == 0) {
					for (int cls = 0; cls < cNum; cls++) {
						data.allNodes[cur_node_id]->fi_parent[cls] = parameter.elnPz[cls];
					}
				}
				else {
					for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {
						int pid = data.allNodes[cur_node_id]->parentsID[p];
						data.allNodes[cur_node_id]->correspondingNeighbour.push_back(pid);
					}
					for (int cls = 0; cls < cNum; cls++) {
						double max = eln(0);
						vector<int> maxNeighbourClass;
						int max_bitCount = 1 << data.allNodes[cur_node_id]->parentsID.size();
						for (int bitCount = 0; bitCount < max_bitCount; bitCount++) { //summation for each parent class label
							vector<int> parentClass;
							double productAccumulator = eln(1);
							int parentClsProd = 1;
							for (int p = 0; p < data.allNodes[cur_node_id]->parentsID.size(); p++) {
								int pid = data.allNodes[cur_node_id]->parentsID[p];
								int parentClsValue = (bitCount >> p) & 1;
								parentClass.push_back(parentClsValue);
								parentClsProd *= parentClsValue;
								productAccumulator = elnproduct(productAccumulator, data.allNodes[pid]->fo[parentClsValue]);  //calculating Product(fo(p)) for all parent of current child except the current node
							}
							productAccumulator = elnproduct(productAccumulator, parameter.elnPz_zpn[cls][parentClsProd]);
							if (max < productAccumulator) {
								max = productAccumulator;
								maxNeighbourClass = parentClass;
							}
							//sumAccumulator = elnsum(sumAccumulator, productAccumulator);
						}
						data.allNodes[cur_node_id]->fi_parent[cls] = max;
						if (cls == 0) {
							for (int t = 0; t < maxNeighbourClass.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassZero.push_back(maxNeighbourClass[t]);
							}
						}
						else {
							for (int t = 0; t < maxNeighbourClass.size(); t++) {
								data.allNodes[cur_node_id]->correspondingNeighbourClassOne.push_back(maxNeighbourClass[t]);
							}
						}
					}
				}

				//calulating fo
				for (int cls = 0; cls < cNum; cls++) { //cls represents class of the current node
					double productAccumulator = eln(1);
					for (int c = 0; c < rightbankchildrenID.size(); c++) {
						int child_id = rightbankchildrenID[c];
						if (child_id == foNode) {
							continue;
						}
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_ChildList[(c * cNum + cls)]); //multiplying with al the child fi except the outgoing child
					}
					productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_parent[cls]);  // multiplying with fi(n)_parent
					productAccumulator = elnproduct(productAccumulator, parameter.elnPzn_xn[(cur_node_id * cNum + cls)]);
					productAccumulator = elnproduct(productAccumulator, -1 * parameter.elnPz[cls]);
					data.allNodes[cur_node_id]->fo[cls] = productAccumulator;
				}

			}

			else {  //message pass n-> parent there is no fi(n)_parent   //computes for root node as well
					//calulating fo
				for (int cls = 0; cls < cNum; cls++) { //cls represents class of the current node
					double productAccumulator = eln(1);
					for (int c = 0; c < rightbankchildrenID.size(); c++) {
						productAccumulator = elnproduct(productAccumulator, data.allNodes[cur_node_id]->fi_ChildList[(c * cNum + cls)]); //multiplying with al the child fi except the outgoing child
					}
					//productAccumulator = elnproduct(productAccumulator, parameter.elnPxn_zn[(cur_node_id*cNum + cls)]);
					productAccumulator = elnproduct(productAccumulator, parameter.elnPzn_xn[(cur_node_id * cNum + cls)]);
					productAccumulator = elnproduct(productAccumulator, -1 * parameter.elnPz[cls]);
					data.allNodes[cur_node_id]->fo[cls] = productAccumulator;
				}
			}

			inferVisited[cur_node_id] = 1;
		}
	}

	updateMapPrediction_right();
}

void cFlood::updateMapPrediction_left() {
	//extracting class
	vector<int> nodeClass(parameter.allPixelSize, -1);
	mappredictions.resize(parameter.orgPixelSize, -1);
	vector<int> visited(parameter.allPixelSize, 0);
	//for left trees 
	std::cout << "Leftbank: " << endl;
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		vector<float> maxcost_accumulator;
		if (!data.hasObservedPixelsLeft[leftOrder] || data.leftbfsOrder[leftOrder].size() < PIXELLIMT) {
			continue;
		}
		float maxCost = MAXCOST;
		std::cout << endl << data.reach_ids[leftOrder] << "->" << data.leftbfsOrder[leftOrder].size() << "->";
		int bfsroot = data.leftbfsRootNodes[leftOrder];
		int orgId;
		if (bfsroot != -1) {
			queue<int> que;
			int nodeCls;
			if (data.allNodes[bfsroot]->fo[0] >= data.allNodes[bfsroot]->fo[1]) {
				nodeCls = 0;
			}
			else {
				nodeCls = 1;
			}
			mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;
			nodeClass[bfsroot] = nodeCls;
			que.push(bfsroot);
			while (!que.empty()) {
				int node = que.front();
				que.pop();
				int nodeCls = nodeClass[node];
				visited[node] = 1;
				for (int c = 0; c < data.allNodes[node]->correspondingNeighbour.size(); c++) {
					int neigh_id = data.allNodes[node]->correspondingNeighbour[c];
					int cClass;
					if (nodeCls == 0) {
						cClass = data.allNodes[node]->correspondingNeighbourClassZero[c];
					}
					else {
						cClass = data.allNodes[node]->correspondingNeighbourClassOne[c];
					}
					nodeClass[neigh_id] = cClass;
					mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls;
					if (!visited[neigh_id]) {
						que.push(neigh_id);
					}

				}
			}
		}

		vector<int> leafnodes;
		//step 1: get list of leave nodes 
		for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
			int nodeId = data.leftbfsOrder[leftOrder][i];
			if (data.allNodes[nodeId]->childrenID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
			if (data.allNodes[nodeId]->parentsID.size() > 1) {
				std::cout << data.allNodes[nodeId]->parentsID.size() << endl;
			}
		}
		//step 2: get chain length
		for (int i = 0; i < leafnodes.size(); i++) {
			int chainLength = 0;
			int nodeId = leafnodes[i];
			int leafnodeId = nodeId;
			int roughness = 1;
			float boundary_cost = -1.0;
			float chainMaxCost = data.allNodes[nodeId]->cost;
			pair<float, float> temp_maxCost_boundaryCost_pair/* = make_pair(chainMaxCost, -1.0)*/;
			pair<int, float> temp_cl_cost_pair /*= make_pair(0, -1.0)*/;
			pair<int, int> temp_chainLength_id/* = make_pair(-1, -1)*/;
			pair<float, int> temp_chainMaxCost_id/* = make_pair(chainMaxCost, -1)*/;
			pair<float, int> temp_cost_boundaryNode_pair /*= make_pair(-1.0,-1)*/;
			pair<int, vector<int>> temp_boundaryNode_leafNodes_pair;
			if (EB == 0) {
				temp_cl_cost_pair = make_pair(0, -1.0);
				temp_chainLength_id = make_pair(-1, -1);
			}
			else if (EB == 1) {
				temp_maxCost_boundaryCost_pair = make_pair(chainMaxCost, -1.0);
				temp_chainMaxCost_id = make_pair(chainMaxCost, -1);
			}
			else if (EB == 2) {
				temp_cost_boundaryNode_pair = make_pair(-1.0, -1);
			}
			//pair<float, float> temp_maxCost_boundaryCost_pair = make_pair(chainMaxCost, -1.0);

			while (data.allNodes[nodeId]->parentsID.size() != 0) {
				if (nodeClass[nodeId] == 1 && boundary_cost < 0.0) {
					bool allchildObserved = true;
					bool observed = true;
					if (BOUNDARY_NODES_OBSERVED == 1) {  // 0 does not exclude any branches// 1 consider pits and tree and unobserved and exclude those branches 
						 //2 excludes branches if the boundary of flood and dry is overlapping with pits layer
						for (int idx = 0; idx < data.allNodes[nodeId]->childrenID.size(); idx++) {
							int childid = data.allNodes[nodeId]->childrenID[idx];
							if (data.allNodes[childid]->isObserved != 1) {
								allchildObserved = false;
								break;
							}
						}

						if (data.allNodes[nodeId]->isObserved == 0) {
							observed = false;
							break;
						}
					}
					if (BOUNDARY_NODES_OBSERVED == 2) {  // uncomment after adding pits and Na identifiers
						//for (int idx = 0; idx < data.allNodes[nodeId]->childrenID.size(); idx++) {
						//	int childid = data.allNodes[nodeId]->childrenID[idx];
						//	if (data.allNodes[childid]->isPits == 1 || data.allNodes[nodeId]->isNa == 1) {
						//		allchildObserved = false;
						//		break;
						//	}
						//}

						//if (data.allNodes[nodeId]->isNa == 1 || data.allNodes[nodeId]->isPits ==1) {
						//	observed = false;
						//}
					}
					if (data.allNodes[nodeId]->roughness == 1 || !observed || !allchildObserved || data.allNodes[nodeId]->isNa == 1) {
						break;
					}
					boundary_cost = data.allNodes[nodeId]->cost;
					if (EB == 0) {
						temp_cl_cost_pair.second = data.allNodes[nodeId]->cost;
						temp_chainLength_id.second = leafnodeId;
					}
					else if (EB == 1) {
						temp_maxCost_boundaryCost_pair.second = data.allNodes[nodeId]->cost;
						temp_chainMaxCost_id.second = leafnodeId;
					}
					else if (EB == 2) {
						temp_cost_boundaryNode_pair.first = data.allNodes[nodeId]->cost;
						temp_cost_boundaryNode_pair.second = nodeId;

						if (data.boundaryNode_leafNodes_Map.find(nodeId) == data.boundaryNode_leafNodes_Map.end()) {
							vector<int> leafs;
							leafs.push_back(leafnodeId);
							data.boundaryNode_leafNodes_Map.insert(make_pair(nodeId, leafs));
						}
						else {
							data.boundaryNode_leafNodes_Map[nodeId].push_back(leafnodeId);
						}
					}
					data.leafNode_boundaryNodes.insert(make_pair(leafnodeId, nodeId));
					roughness = 0;
					if ((EB == 1) || (EB == 2)) {
						break;
					}
				}
				if (EB == 0) {
					chainLength++;
				}
				nodeId = data.allNodes[nodeId]->parentsID[0];
			}
			if (EB == 0) {
				temp_cl_cost_pair.first = chainLength;
				temp_chainLength_id.first = chainLength;
			}
			if (roughness == 0) {
				if (EB == 0) {
					data.chainLength_cost_Pairs.push_back(temp_cl_cost_pair);
					data.chainLength_nodeid_Pairs.push_back(temp_chainLength_id);
				}
				else if (EB == 1) {
					data.maxChainCost_cost_Pairs.push_back(temp_maxCost_boundaryCost_pair);
					data.maxCost_nodeid_Pairs.push_back(temp_chainMaxCost_id);
				}
				else if (EB == 2) {
					data.cost_boundaryNode_Pairs.insert(temp_cost_boundaryNode_pair);
				}
			}
		}
		if (EB == 0) {
			if (data.chainLength_cost_Pairs.size() != data.chainLength_nodeid_Pairs.size()) {
				std::cout << "Length Mismatch Error" << endl;
				exit(0);
			}
		}
		else if (EB == 1) {
			if (data.maxChainCost_cost_Pairs.size() != data.maxCost_nodeid_Pairs.size()) {
				std::cout << "Length Mismatch Error" << endl;
				exit(0);
			}
		}

		if (EB == 0) {
			// using chain length or max cost 
			if (data.chainLength_cost_Pairs.size() != 0) {
				sort(data.chainLength_cost_Pairs.rbegin(), data.chainLength_cost_Pairs.rend());
				sort(data.chainLength_nodeid_Pairs.rbegin(), data.chainLength_nodeid_Pairs.rend());

				//top 20 percent
				//int top = data.chainLength_cost_Pairs.size() * (CUTOFF_PERCENT/100);
				int top = data.chainLength_cost_Pairs.size() * parameter.cutoff_percentage;
				float sum = 0.0;
				for (int j = 0; j < top; j++) {
					sum = sum + data.chainLength_cost_Pairs[j].second;
					if (data.chainLength_nodeid_Pairs[j].second == -1) {
						std::cout << "Issue" << endl;
					}
					//for chain length
					data.boundary_LeafNodes.push_back(data.chainLength_nodeid_Pairs[j].second);

				}
				float avg = -1;
				if (top != 0) {
					avg = sum / top;
				}
				//std::cout << "avg =" << avg << " sum =" << sum << " top = " << top << "data.maxCost_nodeid_Pairs.size = " << data.maxCost_nodeid_Pairs.size() << endl;
				if (avg > parameter.minCost) {
					data.inferedmaxCostLeft[leftOrder] = avg;
				}

				std::cout << data.chainLength_cost_Pairs[0].first << "->" << data.inferedmaxCostLeft[leftOrder] << endl;
			}
		}

		else if (EB == 1) {
			// using max cost branches
			if (data.maxChainCost_cost_Pairs.size() != 0) {

				sort(data.maxChainCost_cost_Pairs.rbegin(), data.maxChainCost_cost_Pairs.rend());
				sort(data.maxCost_nodeid_Pairs.rbegin(), data.maxCost_nodeid_Pairs.rend());


				//top 20 percent
				//int top = data.maxChainCost_cost_Pairs.size() * (CUTOFF_PERCENT/100);
				int top = data.maxChainCost_cost_Pairs.size() * parameter.cutoff_percentage;
				float sum = 0.0;
				for (int j = 0; j < top; j++) {
					sum = sum + data.maxChainCost_cost_Pairs[j].second;
					if (data.maxCost_nodeid_Pairs[j].second == -1) {
						std::cout << "Issue" << endl;
					}
					//for max cost
					data.boundary_LeafNodes.push_back(data.maxCost_nodeid_Pairs[j].second);

				}
				float avg = -1;
				if (top != 0) {
					avg = sum / top;
				}
				//std::cout << "avg =" << avg << " sum =" << sum << " top = " << top << "data.maxCost_nodeid_Pairs.size = " << data.maxCost_nodeid_Pairs.size() << endl;
				if (avg > parameter.minCost) {
					data.inferedmaxCostLeft[leftOrder] = avg;
				}
				if (EB == 1) {
					std::cout << data.maxChainCost_cost_Pairs[0].first << "->" << data.inferedmaxCostLeft[leftOrder] << endl;
				}
				else {
					std::cout << data.chainLength_cost_Pairs[0].first << "->" << data.inferedmaxCostLeft[leftOrder] << endl;
				}

			}
		}

		else if (EB == 2) {
			//using boundary cost
			if (data.cost_boundaryNode_Pairs.size() != 0) {
				if (DEBUG_OUTPUT == 1) {
					vector<float> infered_cost; 
					set<pair<float, int>>::reverse_iterator it;
					for (it = data.cost_boundaryNode_Pairs.rbegin(); it != data.cost_boundaryNode_Pairs.rend(); ++it) {
						pair<float, int> p = *it;
						infered_cost.push_back(p.first);
					}
					//ofstream inferedCosts;
					//inferedCosts.open(CTOutputLocation + "ProfileTables\\"+ to_string(leftOrder) + "_"+ to_string(data.reach_ids[leftOrder])+"_"+parameter.reachId + "_" + "inferedCosts_left.txt");
					////classout.open(CTOutputLocation + CTPrediction);
					//for (int i = 0; i < infered_cost.size(); i++) {
					//	inferedCosts << infered_cost[i] << endl;
					//}
					//inferedCosts.close();
				}

				int top = data.cost_boundaryNode_Pairs.size() * parameter.cutoff_percentage;
				float sum = 0.0;
				set<pair<float, int>>::reverse_iterator it;
				int counter = 0;
				for (it = data.cost_boundaryNode_Pairs.rbegin(); it != data.cost_boundaryNode_Pairs.rend(); ++it) {
					pair<float, int> p = *it;
					int bNode = p.second;
					sum = sum + p.first;
					for (int n = 0; n < data.boundaryNode_leafNodes_Map[bNode].size(); n++) {
						int lNode = data.boundaryNode_leafNodes_Map[bNode][n];
						data.boundary_LeafNodes.push_back(lNode);
					}
					counter++;
					if (counter == top) {
						break;
					}
				}
				float avg = -1.0;
				if (top != 0) {
					avg = sum / top;
				}
				//std::cout << "avg =" << avg << " sum =" << sum << " top = " << top << "data.maxCost_nodeid_Pairs.size = " << data.maxCost_nodeid_Pairs.size() << endl;
				if (avg > parameter.minCost) {
					data.inferedmaxCostLeft[leftOrder] = avg;
				}
				std::cout << "->" << data.inferedmaxCostLeft[leftOrder] << endl;

			}
		}

		data.chainLength_cost_Pairs.clear();
		data.chainLength_cost_Pairs.shrink_to_fit();
		data.chainLength_nodeid_Pairs.clear();
		data.chainLength_nodeid_Pairs.shrink_to_fit();

		data.maxChainCost_cost_Pairs.clear();
		data.maxChainCost_cost_Pairs.shrink_to_fit();
		data.maxCost_nodeid_Pairs.clear();
		data.maxCost_nodeid_Pairs.shrink_to_fit();

		data.cost_boundaryNode_Pairs.clear();
		data.boundaryNode_leafNodes_Map.clear();

	}
}

void cFlood::updateMapPrediction_right() {
	//extracting class
	vector<int> nodeClass(parameter.allPixelSize, -1);
	mappredictions.resize(parameter.orgPixelSize, -1);
	vector<int> visited(parameter.allPixelSize, 0);
	//for right trees 
	std::cout << "Rightbank: " << endl;
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		vector<float> maxcost_accumulator;
		if (!data.hasObservedPixelsRight[rightOrder] || data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
			continue;
		}
		float maxCost = MAXCOST;
		std::cout << endl << data.reach_ids[rightOrder] << "->" << data.rightbfsOrder[rightOrder].size() << "->";
		int bfsroot = data.rightbfsRootNodes[rightOrder];
		int orgId;
		if (bfsroot != -1) {
			queue<int> que;
			int nodeCls;
			if (data.allNodes[bfsroot]->fo[0] >= data.allNodes[bfsroot]->fo[1]) {
				nodeCls = 0;
			}
			else {
				nodeCls = 1;
			}
			mappredictions[data.allNodes[bfsroot]->originalId] = nodeCls;
			nodeClass[bfsroot] = nodeCls;
			que.push(bfsroot);
			while (!que.empty()) {
				int node = que.front();
				que.pop();
				int nodeCls = nodeClass[node];
				visited[node] = 1;
				for (int c = 0; c < data.allNodes[node]->correspondingNeighbour.size(); c++) {
					int neigh_id = data.allNodes[node]->correspondingNeighbour[c];
					int cClass;
					if (nodeCls == 0) {
						cClass = data.allNodes[node]->correspondingNeighbourClassZero[c];
					}
					else {
						cClass = data.allNodes[node]->correspondingNeighbourClassOne[c];
					}
					nodeClass[neigh_id] = cClass;
					mappredictions[data.allNodes[neigh_id]->originalId] = nodeCls;
					if (!visited[neigh_id]) {
						que.push(neigh_id);
					}

				}
			}
		}

		vector<int> leafnodes;
		//step 1: get list of leave nodes 
		for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
			int nodeId = data.rightbfsOrder[rightOrder][i];
			if (data.allNodes[nodeId]->childrenID.size() == 0) {
				leafnodes.push_back(nodeId);
			}
			if (data.allNodes[nodeId]->parentsID.size() > 1) {
				std::cout << data.allNodes[nodeId]->parentsID.size() << endl;
			}
		}
		//step 2: get chain length
		for (int i = 0; i < leafnodes.size(); i++) {
			int chainLength = 0;
			int nodeId = leafnodes[i];
			int leafnodeId = nodeId;
			int roughness = 1;
			float boundary_cost = -1.0;
			float chainMaxCost = data.allNodes[nodeId]->cost;
			pair<float, float> temp_maxCost_boundaryCost_pair/* = make_pair(chainMaxCost, -1.0)*/;
			pair<int, float> temp_cl_cost_pair /*= make_pair(0, -1.0)*/;
			pair<int, int> temp_chainLength_id/* = make_pair(-1, -1)*/;
			pair<float, int> temp_chainMaxCost_id/* = make_pair(chainMaxCost, -1)*/;
			pair<float, int> temp_cost_boundaryNode_pair /*= make_pair(-1.0,-1)*/;
			pair<int, vector<int>> temp_boundaryNode_leafNodes_pair;
			if (EB == 0) {
				temp_cl_cost_pair = make_pair(0, -1.0);
				temp_chainLength_id = make_pair(-1, -1);
			}
			else if (EB == 1) {
				temp_maxCost_boundaryCost_pair = make_pair(chainMaxCost, -1.0);
				temp_chainMaxCost_id = make_pair(chainMaxCost, -1);
			}
			else if (EB == 2) {
				temp_cost_boundaryNode_pair = make_pair(-1.0, -1);
			}
			//pair<float, float> temp_maxCost_boundaryCost_pair = make_pair(chainMaxCost, -1.0);

			while (data.allNodes[nodeId]->parentsID.size() != 0) {
				if (nodeClass[nodeId] == 1 && boundary_cost < 0.0) {
					bool allchildObserved = true;
					bool observed = true;
					if (BOUNDARY_NODES_OBSERVED == 1) {  // 0 does not exclude any branches// 1 consider pits and tree and unobserved and exclude those branches 
						 //2 excludes branches if the boundary of flood and dry is overlapping with pits layer
						for (int idx = 0; idx < data.allNodes[nodeId]->childrenID.size(); idx++) {
							int childid = data.allNodes[nodeId]->childrenID[idx];
							if (data.allNodes[childid]->isObserved != 1) {
								allchildObserved = false;
								break;
							}
						}

						if (data.allNodes[nodeId]->isObserved == 0) {
							observed = false;
							break;
						}
					}
					if (BOUNDARY_NODES_OBSERVED == 2) {  // uncomment after adding pits and Na identifiers
						//for (int idx = 0; idx < data.allNodes[nodeId]->childrenID.size(); idx++) {
						//	int childid = data.allNodes[nodeId]->childrenID[idx];
						//	if (data.allNodes[childid]->isPits == 1 || data.allNodes[nodeId]->isNa == 1) {
						//		allchildObserved = false;
						//		break;
						//	}
						//}

						//if (data.allNodes[nodeId]->isNa == 1 || data.allNodes[nodeId]->isPits ==1) {
						//	observed = false;
						//}
					}
					if (data.allNodes[nodeId]->roughness == 1 || !observed || !allchildObserved || data.allNodes[nodeId]->isNa == 1) {
						break;
					}
					boundary_cost = data.allNodes[nodeId]->cost;
					if (EB == 0) {
						temp_cl_cost_pair.second = data.allNodes[nodeId]->cost;
						temp_chainLength_id.second = leafnodeId;
					}
					else if (EB == 1) {
						temp_maxCost_boundaryCost_pair.second = data.allNodes[nodeId]->cost;
						temp_chainMaxCost_id.second = leafnodeId;
					}
					else if (EB == 2) {
						temp_cost_boundaryNode_pair.first = data.allNodes[nodeId]->cost;
						temp_cost_boundaryNode_pair.second = nodeId;

						if (data.boundaryNode_leafNodes_Map.find(nodeId) == data.boundaryNode_leafNodes_Map.end()) {
							vector<int> leafs;
							leafs.push_back(leafnodeId);
							data.boundaryNode_leafNodes_Map.insert(make_pair(nodeId, leafs));
						}
						else {
							data.boundaryNode_leafNodes_Map[nodeId].push_back(leafnodeId);
						}
					}
					data.leafNode_boundaryNodes.insert(make_pair(leafnodeId, nodeId));
					roughness = 0;
					if ((EB == 1) || (EB == 2)) {
						break;
					}
				}
				if (EB == 0) {
					chainLength++;
				}
				nodeId = data.allNodes[nodeId]->parentsID[0];
			}
			if (EB == 0) {
				temp_cl_cost_pair.first = chainLength;
				temp_chainLength_id.first = chainLength;
			}
			if (roughness == 0) {
				if (EB == 0) {
					data.chainLength_cost_Pairs.push_back(temp_cl_cost_pair);
					data.chainLength_nodeid_Pairs.push_back(temp_chainLength_id);
				}
				else if (EB == 1) {
					data.maxChainCost_cost_Pairs.push_back(temp_maxCost_boundaryCost_pair);
					data.maxCost_nodeid_Pairs.push_back(temp_chainMaxCost_id);
				}
				else if (EB == 2) {
					data.cost_boundaryNode_Pairs.insert(temp_cost_boundaryNode_pair);
				}
			}
		}
		if (EB == 0) {
			if (data.chainLength_cost_Pairs.size() != data.chainLength_nodeid_Pairs.size()) {
				std::cout << "Length Mismatch Error" << endl;
				exit(0);
			}
		}
		else if (EB == 1) {
			if (data.maxChainCost_cost_Pairs.size() != data.maxCost_nodeid_Pairs.size()) {
				std::cout << "Length Mismatch Error" << endl;
				exit(0);
			}
		}

		if (EB == 0) {
			// using chain length or max cost 
			if (data.chainLength_cost_Pairs.size() != 0) {
				sort(data.chainLength_cost_Pairs.rbegin(), data.chainLength_cost_Pairs.rend());
				sort(data.chainLength_nodeid_Pairs.rbegin(), data.chainLength_nodeid_Pairs.rend());

				//top 20 percent
				//int top = data.chainLength_cost_Pairs.size() * (CUTOFF_PERCENT/100);
				int top = data.chainLength_cost_Pairs.size() * parameter.cutoff_percentage;
				float sum = 0.0;
				for (int j = 0; j < top; j++) {
					sum = sum + data.chainLength_cost_Pairs[j].second;
					if (data.chainLength_nodeid_Pairs[j].second == -1) {
						std::cout << "Issue" << endl;
					}
					//for chain length
					data.boundary_LeafNodes.push_back(data.chainLength_nodeid_Pairs[j].second);

				}
				float avg = -1;
				if (top != 0) {
					avg = sum / top;
				}
				//std::cout << "avg =" << avg << " sum =" << sum << " top = " << top << "data.maxCost_nodeid_Pairs.size = " << data.maxCost_nodeid_Pairs.size() << endl;
				if (avg > parameter.minCost) {
					data.inferedmaxCostRight[rightOrder] = avg;
				}

				std::cout << data.chainLength_cost_Pairs[0].first << "->" << data.inferedmaxCostRight[rightOrder] << endl;
			}
		}

		else if (EB == 1) {
			// using max cost branches
			if (data.maxChainCost_cost_Pairs.size() != 0) {

				sort(data.maxChainCost_cost_Pairs.rbegin(), data.maxChainCost_cost_Pairs.rend());
				sort(data.maxCost_nodeid_Pairs.rbegin(), data.maxCost_nodeid_Pairs.rend());


				//top 20 percent
				//int top = data.maxChainCost_cost_Pairs.size() * (CUTOFF_PERCENT/100);
				int top = data.maxChainCost_cost_Pairs.size() * parameter.cutoff_percentage;
				float sum = 0.0;
				for (int j = 0; j < top; j++) {
					sum = sum + data.maxChainCost_cost_Pairs[j].second;
					if (data.maxCost_nodeid_Pairs[j].second == -1) {
						std::cout << "Issue" << endl;
					}
					//for max cost
					data.boundary_LeafNodes.push_back(data.maxCost_nodeid_Pairs[j].second);

				}
				float avg = -1;
				if (top != 0) {
					avg = sum / top;
				}
				//std::cout << "avg =" << avg << " sum =" << sum << " top = " << top << "data.maxCost_nodeid_Pairs.size = " << data.maxCost_nodeid_Pairs.size() << endl;
				if (avg > parameter.minCost) {
					data.inferedmaxCostRight[rightOrder] = avg;
				}
				if (EB == 1) {
					std::cout << data.maxChainCost_cost_Pairs[0].first << "->" << data.inferedmaxCostRight[rightOrder] << endl;
				}
				else {
					std::cout << data.chainLength_cost_Pairs[0].first << "->" << data.inferedmaxCostRight[rightOrder] << endl;
				}

			}
		}

		else if (EB == 2) {
			//using boundary cost
			if (data.cost_boundaryNode_Pairs.size() != 0) {
				if (DEBUG_OUTPUT == 1) {
					vector<float> infered_cost;
					set<pair<float, int>>::reverse_iterator it;
					for (it = data.cost_boundaryNode_Pairs.rbegin(); it != data.cost_boundaryNode_Pairs.rend(); ++it) {
						pair<float, int> p = *it;
						infered_cost.push_back(p.first);
					}
					//ofstream inferedCosts;
					//inferedCosts.open(CTOutputLocation + "ProfileTables\\" + to_string(rightOrder) + "_" + to_string(data.reach_ids[rightOrder]) + "_" + parameter.reachId + "_" + "inferedCosts_right.txt");
					////classout.open(CTOutputLocation + CTPrediction);
					//for (int i = 0; i < infered_cost.size(); i++) {
					//	inferedCosts << infered_cost[i] << endl;
					//}
					//inferedCosts.close();
				}
				int top = data.cost_boundaryNode_Pairs.size() * parameter.cutoff_percentage;
				float sum = 0.0;
				set<pair<float, int>>::reverse_iterator it;
				int counter = 0;
				for (it = data.cost_boundaryNode_Pairs.rbegin(); it != data.cost_boundaryNode_Pairs.rend(); ++it) {
					pair<float, int> p = *it;
					int bNode = p.second;
					sum = sum + p.first;
					for (int n = 0; n < data.boundaryNode_leafNodes_Map[bNode].size(); n++) {
						int lNode = data.boundaryNode_leafNodes_Map[bNode][n];
						data.boundary_LeafNodes.push_back(lNode);
					}
					counter++;
					if (counter == top) {
						break;
					}
				}
				float avg = -1.0;
				if (top != 0) {
					avg = sum / top;
				}
				//std::cout << "avg =" << avg << " sum =" << sum << " top = " << top << "data.maxCost_nodeid_Pairs.size = " << data.maxCost_nodeid_Pairs.size() << endl;
				if (avg > parameter.minCost) {
					data.inferedmaxCostRight[rightOrder] = avg;
				}
				std::cout << "->" << data.inferedmaxCostRight[rightOrder] << endl;

			}
		}

		data.chainLength_cost_Pairs.clear();
		data.chainLength_cost_Pairs.shrink_to_fit();
		data.chainLength_nodeid_Pairs.clear();
		data.chainLength_nodeid_Pairs.shrink_to_fit();

		data.maxChainCost_cost_Pairs.clear();
		data.maxChainCost_cost_Pairs.shrink_to_fit();
		data.maxCost_nodeid_Pairs.clear();
		data.maxCost_nodeid_Pairs.shrink_to_fit();

		data.cost_boundaryNode_Pairs.clear();
		data.boundaryNode_leafNodes_Map.clear();

	}
}

void cFlood::output() {
	auto start = std::chrono::system_clock::now();

	//Note to Saramsha: Save mappredictions to tif
	ofstream classout;
	classout.open(CTOutputLocation + parameter.reachId + "_Prediction.txt");
	//classout.open(CTOutputLocation + parameter.fname + ".txt");
	
	for (int i = 0; i < mappredictions.size(); i++) {
		classout << mappredictions[i] << endl;
	}
	classout.close();
	//Note end
	float** prediction = new float* [parameter.ROW];
	int index = 0;
	for (int row = 0; row < parameter.ROW; row++)
	{
		prediction[row] = new float[parameter.COLUMN];
		for (int col = 0; col < parameter.COLUMN; col++)
		{
			prediction[row][col] = mappredictions[index];
			index++;
		}
	}
	GDALDataset *srcDataset = (GDALDataset*)GDALOpen((CTInputLocation + CTBank).c_str(), GA_ReadOnly);
	double geotransform[6];
	srcDataset->GetGeoTransform(geotransform);
	const OGRSpatialReference* poSRS = srcDataset->GetSpatialRef();
	;
	GeotiffWrite finalTiff((CTOutputLocation + parameter.reachId + "_Prediction.tif").c_str(), parameter.ROW, parameter.COLUMN, 1, geotransform, poSRS);
	finalTiff.write(prediction);
	
	for (int i = 0; i < data.reach_ids.size(); i++) {
		data.avgCost[i] = (data.inferedmaxCostLeft[i] + data.inferedmaxCostRight[i]) / 2.0;
	}
	ofstream profiletable;
	profiletable.open(CTOutputLocation + "ProfileTables\\" + parameter.reachId + "_ProfileTable.csv");
	//profiletable.open(CTOutputLocation + "ProfileTables\\" + parameter.fname + "_ProfileTable.csv");
	profiletable << "SourceId" << "," << "fel" << "," << "Distance Covered" << "," << "Flow Accumulation" << "," << "Left Node Count" << "," << "Left Max Cost" << "," << "Left Infered Cost" << ","
		<< "Observation Indicator" << "," << "Right Node Count" << "," << "Right Max Cost" << "," << "Right Infered Cost" << ","
		<< "Observation Indicator" << "," << "Combined Cost" << "," << "Average Cost" << endl;
	for (int index = 0; index < data.AdjustedReachNodes.size(); index++) {
		profiletable
			<< data.reach_ids[index] << ","
			<< data.allNodes[data.AdjustedReachNodes[index]]->fel << ","
			<< data.distance_covered[index] << ","
			<< data.flow_accumulation[index] << ","
			<< data.leftbfsOrder[index].size() << ","
			<< data.highestCostLeft[index] << ","
			<< data.inferedmaxCostLeft[index] << ","
			<< data.hasObservedPixelsLeft[index] << ","
			<< data.rightbfsOrder[index].size() << ","
			<< data.highestCostRight[index] << ","
			<< data.inferedmaxCostRight[index] << ","
			<< data.hasObservedPixelsRight[index] << ","
			<< data.combinedCost[index] << ","
			<< data.avgCost[index] << endl;
	}
	profiletable.close();

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double>elapsed_seconds = end - start;
	std::cout << "Writing Prediction File took " << elapsed_seconds.count() << "seconds" << endl;
}

int main(int argc, char* argv[]) {
	cFlood flood;
	flood.input(argc, argv);
}


//Testing Functions

void cFlood::getOriginIdBanks() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees 
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[leftOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg;
			}
		}
	}

	//for right trees 
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[rightOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg * 2;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_originIdBanks.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getIds() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees 
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = nodid;
			}
		}
	}

	//for right trees 
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = nodid;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_Ids.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getOrgIds() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees 
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = data.allNodes[nodid]->originalId;
			}
		}
	}

	//for right trees 
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = data.allNodes[nodid]->originalId;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_OrgIds.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getOriginIdLeftBanks() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees 
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[leftOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_originIdLeftBanks.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getOriginIdRightBanks() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for right trees 
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[rightOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_originIdRightBanks.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getRegionNodeCount() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees 
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			int regionSize = data.leftbfsOrder[leftOrder].size();
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = regionSize;
			}
		}
	}

	//for right trees 
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			int regionSize = data.rightbfsOrder[rightOrder].size();
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = regionSize;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_regionSize.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::getLeftRegionNodeCount() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees 
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			int regionSize = data.leftbfsOrder[leftOrder].size();
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = regionSize;
			}
		}
	}

	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_LeftregionSize.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}

void cFlood::getRightRegionNodeCount() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for right trees 
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			int regionSize = data.rightbfsOrder[rightOrder].size();
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = regionSize;
			}
		}
	}
	ofstream classout;
	string reachId = CTPara.substr(9, 2);
	classout.open(CTOutputLocation + reachId + "_RightregionSize.txt");
	for (int i = 0; i < tempMap.size(); i++) {
		classout << tempMap[i] << endl;
	}
	classout.close();
}
void cFlood::sanityChecker() {
	//for left trees 
	//std::cout << "Leftbank: " << endl;
	std::cout << "Starting Sanity Checks..." << endl;
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (!data.hasObservedPixelsLeft[leftOrder] || data.leftbfsOrder[leftOrder].size() < PIXELLIMT) {
			continue;
		}
		// first find list of leaves.
		vector<int> leavesList;
		for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
			int curr = data.leftbfsOrder[leftOrder][i];
			if (data.allNodes[curr]->childrenID.size() == 0) {
				leavesList.push_back(curr);
			}
		}
		for (int i = 0; i < leavesList.size(); i++) {
			int curr = leavesList[i];
			while (data.allNodes[curr]->parentsID.size() != 0) {
				if (mappredictions[data.allNodes[curr]->originalId] == 1) {
					break;
				}
				curr = data.allNodes[curr]->parentsID[0];
			}
			while (data.allNodes[curr]->parentsID.size() != 0) {
				if (mappredictions[data.allNodes[curr]->originalId] != 1) {
					std::cout << i << " !!!Error Elevation Assumption Failed!!!!" << endl;
					break;
				}
				curr = data.allNodes[curr]->parentsID[0];
				//std::cout << "curr= " << curr<<endl;
			}
		}
	}

	//for right trees 
	//std::cout << "Rightbank: " << endl;
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (!data.hasObservedPixelsRight[rightOrder] || data.rightbfsOrder[rightOrder].size() < PIXELLIMT) {
			continue;
		}
		// first find list of leaves.
		vector<int> leavesList;
		for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
			int curr = data.rightbfsOrder[rightOrder][i];
			//std::cout << "curr = " << curr << endl;
			if (data.allNodes[curr]->childrenID.size() == 0) {
				leavesList.push_back(curr);
			}
		}
		for (int i = 0; i < leavesList.size(); i++) {
			int curr = leavesList[i];
			while (data.allNodes[curr]->parentsID.size() != 0) {
				if (mappredictions[data.allNodes[curr]->originalId] == 1) {
					break;
				}
				curr = data.allNodes[curr]->parentsID[0];
			}
			while (data.allNodes[curr]->parentsID.size() != 0) {
				if (mappredictions[data.allNodes[curr]->originalId] != 1) {
					std::cout << i << " !!!Error Elevation Assumption Failed!!!!" << endl;
					break;
				}
				else {
					curr = data.allNodes[curr]->parentsID[0];
				}
			}
		}
	}
	std::cout << "Sanity Test Complete..." << endl;
}

void cFlood::getOriginIdBanks_effectiveBranches() {
	vector<int> tempMap;
	tempMap.resize(parameter.orgPixelSize, -1);
	//for left trees 
	for (int leftOrder = 0; leftOrder < data.leftbfsOrder.size(); leftOrder++) {
		if (data.leftbfsOrder[leftOrder].size() != 0 && data.leftbfsRootNodes[leftOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[leftOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.leftbfsOrder[leftOrder].size(); i++) {
				int nodid = data.leftbfsOrder[leftOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg;
			}
		}
	}

	//for right trees 
	for (int rightOrder = 0; rightOrder < data.rightbfsOrder.size(); rightOrder++) {
		if (data.rightbfsOrder[rightOrder].size() != 0 && data.rightbfsRootNodes[rightOrder] != -1) {
			int reachNode = data.AdjustedReachNodes[rightOrder];
			int reachNodeOrg = data.allNodes[reachNode]->originalId;
			for (int i = 0; i < data.rightbfsOrder[rightOrder].size(); i++) {
				int nodid = data.rightbfsOrder[rightOrder][i];
				tempMap[data.allNodes[nodid]->originalId] = reachNodeOrg * 2;
			}
		}
	}

	for (int i = 0; i < data.boundary_LeafNodes.size(); i++) {
		int nodeId = data.boundary_LeafNodes[i];
		int leafNode = nodeId;
		tempMap[data.allNodes[nodeId]->originalId] = 2;
		nodeId = data.allNodes[nodeId]->parentsID[0];
		//std::cout << "i= " << i << " nodeId = " << nodeId << endl;
		while (data.allNodes[nodeId]->parentsID.size() != 0) {
			tempMap[data.allNodes[nodeId]->originalId] = 1;
			nodeId = data.allNodes[nodeId]->parentsID[0];
		}
		if (data.leafNode_boundaryNodes.find(leafNode) != data.leafNode_boundaryNodes.end()) {
			int bNode = data.leafNode_boundaryNodes[leafNode];
			tempMap[data.allNodes[bNode]->originalId] = 3;
		}
	}

	//ofstream classout;
	//classout.open(CTOutputLocation + parameter.fname + "_Viz2.txt");
	//for (int i = 0; i < tempMap.size(); i++) {
	//	classout << tempMap[i] << endl;
	//}
	//classout.close();
}

