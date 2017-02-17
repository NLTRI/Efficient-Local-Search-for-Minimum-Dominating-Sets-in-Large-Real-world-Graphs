#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <sys/times.h>
#include <unistd.h>
#include <float.h>
#include <iomanip>

using namespace std;

//#define debug_mode

//fixed parameters
const int UNCHANGE_STEP_LIMIT = 100000;
const long long STEP_LIMIT = 10000000000000000; // 10^16
const int BIG_INTEGER = 1000000;
const int TRY_STEP = 100;

//tuned parameters
double lambda;
long long M;

void SAMDS(int cut_off, int vertexNum, vector<int> &DegreeList, vector<vector<int>> &AdjacencyMatrixWithSelfloop, int &best_dominating_set_size, long long &best_solve_step, double &best_cmp_time, double &step_speed)
{
	vector<int> x;
	vector<int> vCoverList;
	vector<double> vCumsumList;
	vector<int> vIsolationList;
	vector<int> neighborList;

	best_dominating_set_size = vertexNum;
	int xSize = 0;
	int xCover = 0;
	double T = -1;
	long long step = 0;
	int counter = 0;

	double randValue;

	int dvIdx;
	int avIdx;
	int nvIdx;

// timing

	tms start, finish;
	int start_time;
	times(&start);
	start_time = start.tms_utime + start.tms_stime;	
	double elap_time = 0;

	/*******************************************************************************/
	// Initialization
	for (int i = 0; i < vertexNum; i++)
	{
		if (DegreeList[i] > 0)
		{
			//randValue = dis(gen);
			randValue = (rand() % BIG_INTEGER) / double(BIG_INTEGER);
			if (randValue >= 0.5)
			{
				x.push_back(1);
				xSize++;
			}				
			else
			{
				x.push_back(0);
			}				
			vIsolationList.push_back(0);
		}
		else
		{
			x.push_back(1);
			xSize++;
			vIsolationList.push_back(1);
		}
		vCoverList.push_back(0);
		vCumsumList.push_back(0);
	}

	// Compute how many times each vertex is covered
	for (int i = 0; i < vertexNum; i++)
	{
		neighborList = AdjacencyMatrixWithSelfloop[i];
		for (int j = 0; j < neighborList.size(); j++)
		{
			nvIdx = neighborList[j];
			vCoverList[i] = vCoverList[i] + x[nvIdx];
		}
		if (vCoverList[i] > 0)
		{
			xCover++;
		}
	}

	/*******************************************************************************/
	// Main loop
	while (true)
	{
		for (int mIdx = 0; mIdx < M; mIdx++)
		{
#ifdef debug_mode
if((step + mIdx) % 100 == 0) cout << "step: " << step + mIdx << endl;
#endif
			//measuring time
			if((step + mIdx) % TRY_STEP == 0)
			{
				times(&finish);
				elap_time = double(finish.tms_utime + finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);
				if(elap_time > cut_off)
				{
					break;
				}			
			}

			double xValue = double(xCover) / vertexNum + 1.0 / (vertexNum * xSize);

			if (xValue >= 1) // Remove vertex
			{
				// Prepare for probabilistic selection
				for (int i = 0; i < vertexNum; i++)
				{
					if (x[i] == 1 && vIsolationList[i] == 0)
					{
						if (i == 0)
							vCumsumList[i] = 1.0 / DegreeList[i];
						else
							vCumsumList[i] = vCumsumList[i - 1] + 1.0 / DegreeList[i];
					}
					else
					{
						if (i == 0)
							vCumsumList[i] = 0;
						else
							vCumsumList[i] = vCumsumList[i - 1];
					}
				}
				for (int i = 0; i < vertexNum; i++)
				{
					vCumsumList[i] = vCumsumList[i] / vCumsumList[vertexNum - 1];
				}
				// Probabilistic selection
				dvIdx = -1;
				//randValue = dis(gen);
				randValue = (rand() % BIG_INTEGER) / double(BIG_INTEGER);				
				for (int i = 0; i < vertexNum; i++)
				{
					if (randValue <= vCumsumList[i])
					{
						dvIdx = i;
						break;
					}
				}
				// Remove and update
				x[dvIdx] = 0;
				xSize--;
				neighborList = AdjacencyMatrixWithSelfloop[dvIdx];
				for (int i = 0; i < neighborList.size(); i++)
				{			
					nvIdx = neighborList[i];
					vCoverList[nvIdx]--;
					if (vCoverList[nvIdx] == 0)
					{
						xCover--;
					}
				}
				// accept or not 
				double yValue = double(xCover) / vertexNum + 1.0 / (vertexNum * xSize);
				double delta = yValue - xValue;
				if (T < 0)
					T = -(abs(delta) / log(0.9));
				if (delta >= 0)
				{
					counter = 0;
					if (xCover == vertexNum && best_dominating_set_size > xSize)
					{
						best_dominating_set_size = xSize;
						times(&finish);
						best_cmp_time = double(finish.tms_utime + finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);
						best_solve_step = step + mIdx;
					}
//
				}
				else
				{
					//randValue = dis(gen);
					randValue = (rand() % BIG_INTEGER) / double(BIG_INTEGER);
					if (randValue < exp(delta / T))
					{
						counter = 0;
						if (xCover == vertexNum && best_dominating_set_size > xSize)
						{
							best_dominating_set_size = xSize;
							times(&finish);
							best_cmp_time = double(finish.tms_utime + finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);
							best_solve_step = step + mIdx;
						}
//
					}
					else
					{
						// Reject and add back
						x[dvIdx] = 1;
						xSize++;
						counter++;
						neighborList = AdjacencyMatrixWithSelfloop[dvIdx];
						for (int i = 0; i < neighborList.size(); i++)
						{
							nvIdx = neighborList[i];
							vCoverList[nvIdx]++;
							if (vCoverList[nvIdx] == 1)
							{
								xCover++;
							}
						}					
					}
				}
			}
			else // Add or swap vertex
			{
				// Prepare for probabilitic selection
				for (int i = 0; i < vertexNum; i++)
				{
					if (x[i] == 0)
					{
						if (i == 0)
							vCumsumList[i] = DegreeList[i];
						else
							vCumsumList[i] = vCumsumList[i - 1] + DegreeList[i];
					}
					else
					{
						if (i == 0)
							vCumsumList[i] = 0;
						else
							vCumsumList[i] = vCumsumList[i - 1];
					}
				}
				for (int i = 0; i < vertexNum; i++)
				{
					vCumsumList[i] = vCumsumList[i] / vCumsumList[vertexNum - 1];
				}
				// Probabilitic selection				
				avIdx = -1;
				//randValue = dis(gen);
				randValue = (rand() % BIG_INTEGER) / double(BIG_INTEGER);
				for (int i = 0; i < vertexNum; i++)
				{
					if (randValue <= vCumsumList[i])
					{
						avIdx = i;
						break;
					}
				}
				// Add
				x[avIdx] = 1;
				xSize++;
				neighborList = AdjacencyMatrixWithSelfloop[avIdx];
				for (int i = 0; i < neighborList.size(); i++)
				{
					nvIdx = neighborList[i];
					vCoverList[nvIdx]++;
					if (vCoverList[nvIdx] == 1)
					{
						xCover++;
					}
				}
				//
				double yValue = double(xCover) / vertexNum + 1.0 / (vertexNum * xSize);
				double delta = yValue - xValue;
				if (T < 0)
					T = -(abs(delta) / log(0.9));
				if (delta > 0)
				{
					counter = 0;
				}
				else
				{
					// reject and remove
					x[avIdx] = 0;
					xSize--;
					neighborList = AdjacencyMatrixWithSelfloop[avIdx];
					for (int i = 0; i < neighborList.size(); i++)
					{
						nvIdx = neighborList[i];
						vCoverList[nvIdx]--;
						if (vCoverList[nvIdx] == 0)
						{
							xCover--;
						}
					}
					// Swap Swap Swap
					// select to Add, with the same cumulative values
					avIdx = -1;
					//randValue = dis(gen);
					randValue = (rand() % BIG_INTEGER) / double(BIG_INTEGER);					
					for (int i = 0; i < vertexNum; i++)
					{
						if (randValue <= vCumsumList[i])
						{
							avIdx = i;
							break;
						}
					}
					// compute new cumulative values for remove
					for (int i = 0; i < vertexNum; i++)
					{
						if (x[i] == 1 && vIsolationList[i] == 0)
						{
							if (i == 0)
								vCumsumList[i] = 1.0 / DegreeList[i];
							else
								vCumsumList[i] = vCumsumList[i - 1] + 1.0 / DegreeList[i];
						}
						else
						{
							if (i == 0)
								vCumsumList[i] = 0;
							else
								vCumsumList[i] = vCumsumList[i - 1];
						}
					}
					for (int i = 0; i < vertexNum; i++)
					{
						vCumsumList[i] = vCumsumList[i] / vCumsumList[vertexNum - 1];
					}					
					dvIdx = -1;
					//randValue = dis(gen);		
					randValue = (rand() % BIG_INTEGER) / double(BIG_INTEGER);			
					for (int i = 0; i < vertexNum; i++)
					{
						if (randValue <= vCumsumList[i])
						{
							dvIdx = i;
							break;
						}
					}
					// Add
					x[avIdx] = 1;
					neighborList = AdjacencyMatrixWithSelfloop[avIdx];
					for (int i = 0; i < neighborList.size(); i++)
					{
						nvIdx = neighborList[i];
						vCoverList[nvIdx]++;
						if (vCoverList[nvIdx] == 1)
						{
							xCover++;
						}
					}
					// Remove
					x[dvIdx] = 0;
					neighborList = AdjacencyMatrixWithSelfloop[dvIdx];
					for (int i = 0; i < neighborList.size(); i++)
					{
						nvIdx = neighborList[i];
						vCoverList[nvIdx]--;
						if (vCoverList[nvIdx] == 0)
						{
							xCover--;
						}
					}					
					double yValue = double(xCover) / vertexNum + 1.0 / (vertexNum * xSize);
					double delta = yValue - xValue;
					if (T < 0)
						T = -(abs(delta) / log(0.9));
					if (delta >= 0)
					{
						counter = 0;
						if (xCover == vertexNum && best_dominating_set_size > xSize)
						{
							best_dominating_set_size = xSize;
							times(&finish);
							best_cmp_time = double(finish.tms_utime + finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);
							best_solve_step = step + mIdx;
						}
					}
					else
					{
						//randValue = dis(gen);
						randValue = (rand() % BIG_INTEGER) / double(BIG_INTEGER);
						if (randValue < exp(delta / T))
						{
							counter = 0;
							if (xCover == vertexNum && best_dominating_set_size > xSize)
							{
								best_dominating_set_size = xSize;
								times(&finish);
								best_cmp_time = double(finish.tms_utime + finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);
								best_solve_step = step + mIdx;
							}
						}
						else
						{
							// Restore
							x[avIdx] = 0;
							neighborList = AdjacencyMatrixWithSelfloop[avIdx];
							for (int i = 0; i < neighborList.size(); i++)
							{
								nvIdx = neighborList[i];
								vCoverList[nvIdx]--;
								if (vCoverList[nvIdx] == 0)
								{
									xCover--;
								}
							}
							x[dvIdx] = 1;	
							neighborList = AdjacencyMatrixWithSelfloop[dvIdx];
							for (int i = 0; i < neighborList.size(); i++)
							{
								nvIdx = neighborList[i];
								vCoverList[nvIdx]++;
								if (vCoverList[nvIdx] == 1)
								{
									xCover++;
								}
							}
							counter++;
						}
					}
				}

			}
		}

		step = step + M;

		if (step >= STEP_LIMIT || counter >= UNCHANGE_STEP_LIMIT || elap_time > cut_off)
			break;
		else
			T = T * lambda;
	}

	times(&finish);
	elap_time = double(finish.tms_utime + finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);

	step_speed = (long double)(step) / 1000.0 / elap_time;
}

int main(int argc, char* argv[])
{
	char filename[1024];
	sscanf(argv[1], "%s", filename);

	int seed;
	sscanf(argv[2], "%d", &seed);
	srand(seed);

	sscanf(argv[3], "%lf", &lambda);
	sscanf(argv[4], "%lld", &M);

	int cut_off;//number of seconds
	sscanf(argv[5], "%d", &cut_off);

	int best_dominating_set_size;
	long long best_solve_step = 0;
	double best_cmp_time = 0;
	double step_speed = DBL_MAX;
	/****************************************************************************************/
/*
	ifstream infile;
	infile.open(filename);
*/
	ifstream infile(filename);
	if (!infile)
	{
		return 1;
	}

	char line[1024];
	infile.getline(line, 1024);
	while (line[0] != 'p')
		infile.getline(line, 1024);
	int v_num, e_num;
	char tmpStr1[1024], tmpStr2[1024];
	sscanf(line, "%s %s %d %d", tmpStr1, tmpStr2, &v_num, &e_num);

	vector<vector<int>> AdjacencyMatrixWithSelfloop;
	vector<int> DegreeList;
	for (int i = 0; i < v_num; i++)
	{
		vector<int> neighborList;
		neighborList.push_back(i);
		AdjacencyMatrixWithSelfloop.push_back(neighborList);
		DegreeList.push_back(0);
	}

	char tmpChar;
	int v1, v2;	
	for (int i = 0; i < e_num; i++)
	{
		infile >> tmpChar >> v1 >> v2;
		AdjacencyMatrixWithSelfloop[v1 - 1].push_back(v2 - 1);
		AdjacencyMatrixWithSelfloop[v2 - 1].push_back(v1 - 1);
		DegreeList[v1 - 1] = DegreeList[v1 - 1] + 1;
		DegreeList[v2 - 1] = DegreeList[v2 - 1] + 1;
	}

	infile.close();

#ifdef debug_mode
cout << "the problem instance has been read successfully" << endl;
#endif
	
	/****************************************************************************************/
	SAMDS(cut_off, v_num, DegreeList, AdjacencyMatrixWithSelfloop, best_dominating_set_size, best_solve_step, best_cmp_time, step_speed);

	cout << "o " << best_dominating_set_size << endl;
	cout << "c solveTime " << best_cmp_time << endl;
	cout << "c solveStep " << best_solve_step << endl;
	cout << "c stepSpeed(/ms) " << fixed << setprecision(2) << step_speed << endl;
	return 0;
}
