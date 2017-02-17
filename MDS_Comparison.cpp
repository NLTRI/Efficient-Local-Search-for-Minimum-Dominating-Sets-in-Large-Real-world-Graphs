#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>  
#include <cmath>
#include <omp.h>
#include <chrono>
#include <ctime>

using namespace std;
using get_time = chrono::steady_clock;

vector<int> sortDegreeVectorAscendingIndex(const vector<int> &degreeVector)
{
	// initialize original index locations
	vector<int> degreeAscendIndexVector(degreeVector.size());
	iota(degreeAscendIndexVector.begin(), degreeAscendIndexVector.end(), 0);

	// sort indexes based on comparing values in degreeVector
	sort(degreeAscendIndexVector.begin(), degreeAscendIndexVector.end(),
		[&degreeVector](int i1, int i2) {return degreeVector[i1] < degreeVector[i2]; });

	return degreeAscendIndexVector;
}

vector<int> sortDegreeVectorDescendingIndex(const vector<int> &degreeVector)
{
	// initialize original index locations
	vector<int> degreeDescendIndexVector(degreeVector.size());
	iota(degreeDescendIndexVector.begin(), degreeDescendIndexVector.end(), 0);

	// sort indexes based on comparing values in degreeVector
	sort(degreeDescendIndexVector.begin(), degreeDescendIndexVector.end(),
		[&degreeVector](int i1, int i2) {return degreeVector[i1] > degreeVector[i2]; });

	return degreeDescendIndexVector;
}

int Algorithm_1(const vector<vector<int>* > &coverMatrix, const vector<int> &degreeVector)
{
	vector<bool> mdsIndicator;
	vector<bool> coverIndicator;
	vector<int> localCoverageGain;
	vector<int> *thisCoverIndex, *thisNeighborCoverIndex;
	int coveredNum = 0;
	vector<bool> maxValueVertexIndicator;
	vector<int> maxValueVertexVector;
	int maxValue, maxIndex;
	int thisNeighborIndex, jIdx;

	int vertexNum = coverMatrix.size();

	for (int i = 0; i < vertexNum; i++)
	{
		if (degreeVector[i] == 0)
		{
			mdsIndicator.push_back(true);
			coverIndicator.push_back(true);
			coveredNum++;
		}
		else
		{
			mdsIndicator.push_back(false);
			coverIndicator.push_back(false);
		}
		localCoverageGain.push_back(0);
		maxValueVertexIndicator.push_back(false);
	}

	for (int i = 0; i < vertexNum; i++)
	{
		if (!mdsIndicator[i])
		{
			thisCoverIndex = coverMatrix[i];
			for (int j = 0; j < thisCoverIndex->size(); j++)
			{
				if (!coverIndicator[thisCoverIndex->at(j)])
				{
					localCoverageGain[i] = localCoverageGain[i] + 1;
				}
			}
		}
	}

	while (coveredNum < vertexNum)
	{
		maxValue = 0;
		
		for (int i = 0; i < vertexNum; i++)
		{
			if (!mdsIndicator[i])
			{				
				if (maxValue < localCoverageGain[i])
				{
					maxValue = localCoverageGain[i];
				}
			}
			maxValueVertexIndicator[i] = false;
		}
		maxValueVertexVector.clear();
		for (int i = 0; i < vertexNum; i++)
		{
			if (!mdsIndicator[i] && localCoverageGain[i] == maxValue)
			{
				maxValueVertexIndicator[i] = true;
				maxValueVertexVector.push_back(i);
			}
		}

		for (int mvIdx = 0; mvIdx < maxValueVertexVector.size(); mvIdx++)
		{
			if (coveredNum == vertexNum)
			{
				break;
			}
			maxIndex = maxValueVertexVector[mvIdx];
			if (maxValueVertexIndicator[maxIndex])
			{
				mdsIndicator[maxIndex] = true;
				thisCoverIndex = coverMatrix[maxIndex];
				for (int i = 0; i < thisCoverIndex->size(); i++)
				{
					thisNeighborIndex = thisCoverIndex->at(i);
					if (!coverIndicator[thisNeighborIndex])
					{
						coverIndicator[thisNeighborIndex] = true;
						coveredNum++;
						thisNeighborCoverIndex = coverMatrix[thisNeighborIndex];
						for (int j = 0; j < thisNeighborCoverIndex->size(); j++)
						{
							jIdx = thisNeighborCoverIndex->at(j);
							if (!mdsIndicator[jIdx])
							{
								localCoverageGain[jIdx] = localCoverageGain[jIdx] - 1;
								maxValueVertexIndicator[jIdx] = false;
							}
						}
					}
				}
			}
		}
	}

	int mdsSize = 0;
	for (int i = 0; i < vertexNum; i++)
	{
		if (mdsIndicator[i])
		{
			mdsSize++;
		}
	}

	return mdsSize;
}

int Algorithm_2(const vector<vector<int>* > &coverMatrix, const vector<int> &degreeVector)
{
	vector<bool> mdsIndicator;
	vector<bool> coverIndicator;
	vector<int> localCoverageGain;
	vector<int> *thisCoverIndex, *thisNeighborCoverIndex;
	int coveredNum = 0;
	vector<bool> maxValueVertexIndicator;
	vector<int> maxValueVertexVector;
	int maxValue, maxIndex;
	int thisNeighborIndex, jIdx;

	int vertexNum = coverMatrix.size();

	for (int i = 0; i < vertexNum; i++)
	{
		if (degreeVector[i] == 0)
		{
			mdsIndicator.push_back(true);
			coverIndicator.push_back(true);
			coveredNum++;
		}
		else
		{
			mdsIndicator.push_back(false);
			coverIndicator.push_back(false);
		}
		localCoverageGain.push_back(0);
		maxValueVertexIndicator.push_back(false);
	}
	for (int i = 0; i < vertexNum; i++)
	{
		if (degreeVector[i] == 1)
		{
			thisNeighborIndex = coverMatrix[i]->at(1);
			mdsIndicator[thisNeighborIndex] = true;
			mdsIndicator[i] = false;
			thisCoverIndex = coverMatrix[thisNeighborIndex];
			for (int j = 0; j < thisCoverIndex->size(); j++)
			{
				if (!coverIndicator[thisCoverIndex->at(j)])
				{
					coverIndicator[thisCoverIndex->at(j)] = true;
					coveredNum++;
				}
			}
		}
	}

	for (int i = 0; i < vertexNum; i++)
	{
		if (!mdsIndicator[i])
		{
			thisCoverIndex = coverMatrix[i];
			for (int j = 0; j < thisCoverIndex->size(); j++)
			{
				if (!coverIndicator[thisCoverIndex->at(j)])
				{
					localCoverageGain[i] = localCoverageGain[i] + 1;
				}
			}
		}
	}

	while (coveredNum < vertexNum)
	{
		maxValue = 0;

		for (int i = 0; i < vertexNum; i++)
		{
			if (!mdsIndicator[i])
			{
				if (maxValue < localCoverageGain[i])
				{
					maxValue = localCoverageGain[i];
				}
			}
			maxValueVertexIndicator[i] = false;
		}
		maxValueVertexVector.clear();
		for (int i = 0; i < vertexNum; i++)
		{
			if (!mdsIndicator[i] && localCoverageGain[i] == maxValue)
			{
				maxValueVertexIndicator[i] = true;
				maxValueVertexVector.push_back(i);
			}
		}

		for (int mvIdx = 0; mvIdx < maxValueVertexVector.size(); mvIdx++)
		{
			if (coveredNum == vertexNum)
			{
				break;
			}
			maxIndex = maxValueVertexVector[mvIdx];
			if (maxValueVertexIndicator[maxIndex])
			{
				mdsIndicator[maxIndex] = true;
				thisCoverIndex = coverMatrix[maxIndex];
				for (int i = 0; i < thisCoverIndex->size(); i++)
				{
					thisNeighborIndex = thisCoverIndex->at(i);
					if (!coverIndicator[thisNeighborIndex])
					{
						coverIndicator[thisNeighborIndex] = true;
						coveredNum++;
						thisNeighborCoverIndex = coverMatrix[thisNeighborIndex];
						for (int j = 0; j < thisNeighborCoverIndex->size(); j++)
						{
							jIdx = thisNeighborCoverIndex->at(j);
							if (!mdsIndicator[jIdx])
							{
								localCoverageGain[jIdx] = localCoverageGain[jIdx] - 1;
								maxValueVertexIndicator[jIdx] = false;
							}
						}
					}
				}
			}
		}
	}

	int mdsSize = 0;
	for (int i = 0; i < vertexNum; i++)
	{
		if (mdsIndicator[i])
		{
			mdsSize++;
		}
	}

	return mdsSize;
}

int Algorithm_3(const vector<vector<int>* > &coverMatrix, const vector<int> &degreeVector)
{
	vector<bool> mdsIndicator;
	vector<bool> coverIndicator;
	vector<int> localCoverageGain;
	vector<int> *thisCoverIndex, *thisNeighborCoverIndex;
	int coveredNum = 0;
	vector<bool> maxValueVertexIndicator;
	vector<int> maxValueVertexVector;
	int maxValue, maxIndex;
	int thisNeighborIndex, jIdx;

	int vertexNum = coverMatrix.size();

	for (int i = 0; i < vertexNum; i++)
	{
		if (degreeVector[i] == 0)
		{
			mdsIndicator.push_back(true);
			coverIndicator.push_back(true);
			coveredNum++;
		}
		else
		{
			mdsIndicator.push_back(false);
			coverIndicator.push_back(false);
		}
		localCoverageGain.push_back(0);
		maxValueVertexIndicator.push_back(false);
	}

	for (int i = 0; i < vertexNum; i++)
	{
		if (!mdsIndicator[i])
		{
			thisCoverIndex = coverMatrix[i];
			for (int j = 0; j < thisCoverIndex->size(); j++)
			{
				if (!coverIndicator[thisCoverIndex->at(j)])
				{
					localCoverageGain[i] = localCoverageGain[i] + 1;
				}
			}
		}
	}

	while (coveredNum < vertexNum)
	{
		maxValue = 0;

		for (int i = 0; i < vertexNum; i++)
		{
			if (!coverIndicator[i])
			{
				if (maxValue < localCoverageGain[i])
				{
					maxValue = localCoverageGain[i];
				}
			}
			maxValueVertexIndicator[i] = false;
		}
		maxValueVertexVector.clear();
		for (int i = 0; i < vertexNum; i++)
		{
			if (!coverIndicator[i] && localCoverageGain[i] == maxValue)
			{
				maxValueVertexIndicator[i] = true;
				maxValueVertexVector.push_back(i);
			}
		}

		for (int mvIdx = 0; mvIdx < maxValueVertexVector.size(); mvIdx++)
		{
			if (coveredNum == vertexNum)
			{
				break;
			}
			maxIndex = maxValueVertexVector[mvIdx];
			if (maxValueVertexIndicator[maxIndex])
			{
				mdsIndicator[maxIndex] = true;
				thisCoverIndex = coverMatrix[maxIndex];
				for (int i = 0; i < thisCoverIndex->size(); i++)
				{
					thisNeighborIndex = thisCoverIndex->at(i);
					if (!coverIndicator[thisNeighborIndex])
					{
						coverIndicator[thisNeighborIndex] = true;
						coveredNum++;
						thisNeighborCoverIndex = coverMatrix[thisNeighborIndex];
						for (int j = 0; j < thisNeighborCoverIndex->size(); j++)
						{
							jIdx = thisNeighborCoverIndex->at(j);
							if (!coverIndicator[jIdx])
							{
								localCoverageGain[jIdx] = localCoverageGain[jIdx] - 1;
								maxValueVertexIndicator[jIdx] = false;
							}
						}
					}
				}
			}
		}
	}

	int mdsSize = 0;
	for (int i = 0; i < vertexNum; i++)
	{
		if (mdsIndicator[i])
		{
			mdsSize++;
		}
	}

	return mdsSize;
}

int Algorithm_4(const vector<vector<int>* > &coverMatrix, const vector<int> &degreeVector)
{
	vector<bool> mdsIndicator;
	vector<bool> coverIndicator;
	vector<int> localCoverageGain;
	vector<int> *thisCoverIndex, *thisNeighborCoverIndex;
	int coveredNum = 0;
	vector<bool> maxValueVertexIndicator;
	vector<int> maxValueVertexVector;
	int maxValue, maxIndex;
	int thisNeighborIndex, jIdx;

	int vertexNum = coverMatrix.size();

	for (int i = 0; i < vertexNum; i++)
	{
		if (degreeVector[i] == 0)
		{
			mdsIndicator.push_back(true);
			coverIndicator.push_back(true);
			coveredNum++;
		}
		else
		{
			mdsIndicator.push_back(false);
			coverIndicator.push_back(false);
		}
		localCoverageGain.push_back(0);
		maxValueVertexIndicator.push_back(false);
	}
	for (int i = 0; i < vertexNum; i++)
	{
		if (degreeVector[i] == 1)
		{
			thisNeighborIndex = coverMatrix[i]->at(1);
			mdsIndicator[thisNeighborIndex] = true;
			mdsIndicator[i] = false;
			thisCoverIndex = coverMatrix[thisNeighborIndex];
			for (int j = 0; j < thisCoverIndex->size(); j++)
			{
				if (!coverIndicator[thisCoverIndex->at(j)])
				{
					coverIndicator[thisCoverIndex->at(j)] = true;
					coveredNum++;
				}
			}
		}
	}

	for (int i = 0; i < vertexNum; i++)
	{
		if (!mdsIndicator[i])
		{
			thisCoverIndex = coverMatrix[i];
			for (int j = 0; j < thisCoverIndex->size(); j++)
			{
				if (!coverIndicator[thisCoverIndex->at(j)])
				{
					localCoverageGain[i] = localCoverageGain[i] + 1;
				}
			}
		}
	}

	while (coveredNum < vertexNum)
	{
		maxValue = 0;

		for (int i = 0; i < vertexNum; i++)
		{
			if (!coverIndicator[i])
			{
				if (maxValue < localCoverageGain[i])
				{
					maxValue = localCoverageGain[i];
				}
			}
			maxValueVertexIndicator[i] = false;
		}
		maxValueVertexVector.clear();
		for (int i = 0; i < vertexNum; i++)
		{
			if (!coverIndicator[i] && localCoverageGain[i] == maxValue)
			{
				maxValueVertexIndicator[i] = true;
				maxValueVertexVector.push_back(i);
			}
		}

		for (int mvIdx = 0; mvIdx < maxValueVertexVector.size(); mvIdx++)
		{
			if (coveredNum == vertexNum)
			{
				break;
			}
			maxIndex = maxValueVertexVector[mvIdx];
			if (maxValueVertexIndicator[maxIndex])
			{
				mdsIndicator[maxIndex] = true;
				thisCoverIndex = coverMatrix[maxIndex];
				for (int i = 0; i < thisCoverIndex->size(); i++)
				{
					thisNeighborIndex = thisCoverIndex->at(i);
					if (!coverIndicator[thisNeighborIndex])
					{
						coverIndicator[thisNeighborIndex] = true;
						coveredNum++;
						thisNeighborCoverIndex = coverMatrix[thisNeighborIndex];
						for (int j = 0; j < thisNeighborCoverIndex->size(); j++)
						{
							jIdx = thisNeighborCoverIndex->at(j);
							if (!coverIndicator[jIdx])
							{
								localCoverageGain[jIdx] = localCoverageGain[jIdx] - 1;
								maxValueVertexIndicator[jIdx] = false;
							}
						}
					}
				}
			}
		}
	}

	int mdsSize = 0;
	for (int i = 0; i < vertexNum; i++)
	{
		if (mdsIndicator[i])
		{
			mdsSize++;
		}
	}

	return mdsSize;
}

int Algorithm_5(const vector<vector<int>* > &coverMatrix, const vector<int> &degreeVector, vector<int> &degreeDescendIndexVector)
{
	vector<bool> mdsIndicator;
	vector<bool> coverIndicator;
	vector<int> *thisCoverIndex;
	int coveredNum = 0;
	bool isUseful = false;

	int vertexNum = coverMatrix.size();

	for (int i = 0; i < vertexNum; i++)
	{
		if (degreeVector[i] == 0)
		{
			mdsIndicator.push_back(true);
			coverIndicator.push_back(true);
			coveredNum++;
		}
		else
		{
			mdsIndicator.push_back(false);
			coverIndicator.push_back(false);
		}
	}
	
	for (int i = 0; i < vertexNum; i++)
	{
		if (coveredNum == vertexNum)
		{
			break;
		}

		int thisVertex = degreeDescendIndexVector[i];
		if (!mdsIndicator[thisVertex])
		{
			isUseful = false;
			thisCoverIndex = coverMatrix[thisVertex];

			for (int j = 0; j < thisCoverIndex->size(); j++)
			{
				if (!coverIndicator[thisCoverIndex->at(j)])
				{
					coverIndicator[thisCoverIndex->at(j)] = true;
					coveredNum++;
					isUseful = true;
				}
			}
			if (isUseful)
			{
				mdsIndicator[thisVertex] = true;
			}
		}
	}

	int mdsSize = 0;
	for (int i = 0; i < vertexNum; i++)
	{
		if (mdsIndicator[i])
		{
			mdsSize++;
		}
	}

	return mdsSize;
}

int main(int argc, char *argv[])
{
	mt19937 g(static_cast<uint32_t>(time(0)));

	ifstream inFile;
	ofstream outFile;
	stringstream ss;
	string line, token;

	// read all data file path	
	inFile.open(argv[1]);
	if (!inFile)
	{
		return 0;
	}
	vector<string> filePathVector;
	while (getline(inFile, line))
	{
		filePathVector.push_back(line);
	}
	inFile.close();
	int fileIndexStart = atoi(argv[2]); // zero-based
	int fileIndexEnd = atoi(argv[3]);
	int fileNum = fileIndexEnd - fileIndexStart + 1;

	int trialNum = atoi(argv[4]);
	int threadNum = atoi(argv[5]);

	// create containers of results
	int algrorithmNum = 5;
	vector<vector<vector<int> > > result;
	for (int i = 0; i < fileNum; i++)
	{
		vector<vector<int> > thisVectorVector;
		for (int j = 0; j < algrorithmNum; j++)
		{
			vector<int> thisVector;
			for (int k = 0; k < trialNum; k++)
			{
				thisVector.push_back(0);
			}
			thisVectorVector.push_back(thisVector);
		}
		result.push_back(thisVectorVector);
	}

	/*************************************************************************************************************/
	// process each graph
	for (int fileIdx = 0; fileIdx < fileNum; fileIdx++)
	{	
		cout << fileIdx << endl;
		auto start = get_time::now();

		// read data
		string filePath = filePathVector[fileIndexStart + fileIdx];
		if (fileIndexStart + fileIdx < filePathVector.size() - 1)
		{
			filePath.pop_back();
		}
		inFile.open(filePath);
		if (!inFile)
		{
			return 0;
		}

		getline(inFile, line);
		ss.clear();
		ss.str(line);
		getline(ss, token, ' ');
		getline(ss, token, ' ');
		getline(ss, token, ' ');
		int vertexNum = stoi(token);
		getline(ss, token, ' ');
		int edgeNum = stoi(token);

		vector<vector<int> > coverMatrix;
		for (int i = 0; i < vertexNum; i++)
		{
			vector<int> coverVector;
			coverVector.push_back(i);
			coverMatrix.push_back(coverVector);
		}

		int v1, v2;
		while (getline(inFile, line))
		{
			ss.clear();
			ss.str(line);
			getline(ss, token, ' ');
			getline(ss, token, ' ');
			v1 = stoi(token);
			getline(ss, token, ' ');
			v2 = stoi(token);
			coverMatrix[v1 - 1].push_back(v2 - 1);
			coverMatrix[v2 - 1].push_back(v1 - 1);
		}
		inFile.close();


		#pragma omp parallel for num_threads(threadNum)
		for (int trialIdx = 0; trialIdx < trialNum; trialIdx++)
		{
			// random shuffle coverMatrix
			vector<vector<int>* > thisCoverMatrix;
			vector<int> indexVector;
			for (int i = 0; i < vertexNum; i++)
			{
				indexVector.push_back(i);
				vector<int> *coverVector = new vector<int>;
				thisCoverMatrix.push_back(coverVector);
			}

			/*************************************************************************************/

			shuffle(indexVector.begin(), indexVector.end(), g);

			/*************************************************************************************/

			for (int i = 0; i < vertexNum; i++)
			{
				vector<int> coverVector = coverMatrix[i];
				for (int j = 0; j < coverVector.size(); j++)
				{
					thisCoverMatrix[indexVector[i]]->push_back(indexVector[coverVector[j]]);
				}
			}

			vector<int> degreeVector;
			for (int i = 0; i < vertexNum; i++)
			{
				degreeVector.push_back(thisCoverMatrix[i]->size() - 1);
			}
			vector<int> degreeDescendIndexVector = sortDegreeVectorDescendingIndex(degreeVector);
			
			result[fileIdx][0][trialIdx] = Algorithm_1(thisCoverMatrix, degreeVector);
			result[fileIdx][1][trialIdx] = Algorithm_2(thisCoverMatrix, degreeVector);
			result[fileIdx][2][trialIdx] = Algorithm_3(thisCoverMatrix, degreeVector);
			result[fileIdx][3][trialIdx] = Algorithm_4(thisCoverMatrix, degreeVector);
			result[fileIdx][4][trialIdx] = Algorithm_5(thisCoverMatrix, degreeVector, degreeDescendIndexVector);

			for (int i = 0; i < vertexNum; i++)
			{
				delete thisCoverMatrix[i];
			}
		}

		auto end = get_time::now();
		auto diff = end - start;
		cout << chrono::duration_cast<chrono::milliseconds>(diff).count() / 1000.0 << endl;

		vector<double> avgResult;
		for (int j = 0; j < algrorithmNum; j++)
		{
			avgResult.push_back(0);
		}
		for (int j = 0; j < algrorithmNum; j++)
		{
			for (int k = 0; k < trialNum; k++)
			{
				avgResult[j] = avgResult[j] + result[fileIdx][j][k];
			}
			avgResult[j] = avgResult[j] / trialNum;
			cout << avgResult[j] << " ";
		}
		cout << endl;
	}

	/*************************************************************************************************************/
	// output results
	vector<vector<double> > avgResult;
	for (int i = 0; i < fileNum; i++)
	{
		vector<double> thisVector;
		for (int j = 0; j < algrorithmNum; j++)
		{
			thisVector.push_back(0);
		}
		avgResult.push_back(thisVector);
	}
	for (int i = 0; i < fileNum; i++)
	{
		for (int j = 0; j < algrorithmNum; j++)
		{
			for (int k = 0; k < trialNum; k++)
			{
				avgResult[i][j] = avgResult[i][j] + result[i][j][k];
			}
			avgResult[i][j] = avgResult[i][j] / trialNum;
		}
	}
	string resultFileName = "result_Comparison_" + to_string(fileIndexStart) + "_" + to_string(fileIndexEnd) + ".txt";
	outFile.open(resultFileName);
	for (int i = 0; i < avgResult.size(); i++)
	{
		vector<double> thisVector = avgResult[i];
		for (int j = 0; j < thisVector.size(); j++)
		{
			outFile << thisVector[j] << " ";
		}
		outFile << "\n";
	}
	outFile.close();
}