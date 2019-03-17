#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <mpi.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <filesystem>


#define TIME

using namespace std;


// Global Proc Info:
int id, p;
int proc_rows, proc_cols;		// Processor domain dimensions (proc_rows x proc_cols)

// Indv Proc Info:
int id_row, id_col;				// Position of the processor in the processor domain
int id_xsize, id_ysize;			// Number of cells in the processor

// Simulation Settings:
int imax = 100, jmax = 100;		// Cell domain dimensions (imax x jmax)
int max_steps = 20;				// Time steps
int periodic = 1;				// Periodic = 1,   Non-Periodic = 0


double wTime()
{
	return (double)clock() / CLOCKS_PER_SEC;	// Return wall time
}


void find_dimensions(int p, int &rows, int &columns)
{
	int min_gap = p;
	int top = sqrt(p) + 1;									// Most efficient for a square domain.

	for (int i = 1; i <= top; i++)							// for(int i = 1; i <= p / 2; i++)				// Can do more efficiently if sqr rooted.
	{
		if (p%i == 0)										// i is a factor of p
		{
			int gap = abs(p / i - i);						// gap between the factors

			if (gap < min_gap)
			{
				min_gap = gap;
				rows = i;
				columns = p / i;
			}
		}
	}

#ifndef TIME
	if (id == 0)
		cout << "Divide " << p << " into " << rows << " by " << columns << " grid" << endl;
#endif // !TIME

}


void id_to_index(int id, int &id_row, int &id_cols)
{
	id_cols = id % proc_cols;
	id_row = id / proc_cols;
}


int id_from_index(int id_row, int id_col)
{
	if (id_row >= proc_rows || id_row < 0)
		return -1;
	if (id_col >= proc_cols || id_col < 0)
		return -1;

	return id_row * proc_cols + id_col;
}


void num_cells(int &ijmax, int &id_pos, int &id_size, int &proc)
{
	int remaining = ijmax;
	for (int i = 0; i < proc; i++)
	{
		int allocated = remaining / (proc - i);

		if (i <= id_pos)
		{
			id_size = allocated;
			remaining -= allocated;
		}
		else break;
	}
}


int num_neighbours(bool**grid, int ii, int jj)	//ii, jj centre coordinate
{
	int ix, jx;
	int cnt = 0;
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
			if (i != 0 || j != 0)	// if not in the centre coordinate
			{
				ix = i + ii;
				jx = j + jj; 
				if (grid[ix][jx]) cnt++;
			}
	return cnt;
}


void grid_to_file(int it, string main_dir, bool **grid, int id_xsize, int id_ysize)
{
	if (!filesystem::exists(main_dir + "/" + to_string(it)))
	{
		filesystem::create_directory(main_dir + "/" + to_string(it));
	}

	stringstream fname;
	fstream f1;

	fname << main_dir + "/" << it << "/" << id_row << "x" << id_col << ".csv";
		
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 1; i < id_ysize + 1; i++)
	{
		for (int j = 1; j < id_xsize + 1; j++)
			f1 << grid[i][j] << ",";
		f1 << endl;
	}
	f1.close();
}


void do_iteration(bool **grid, bool **new_grid, int id_xsize, int id_ysize)
{
	for (int i = 0; i < id_ysize + 2; i++)
	{
		for (int j = 0; j < id_xsize + 2; j++)
		{
			new_grid[i][j] = grid[i][j];

			if (i >= 1 && i <= id_ysize && j >= 1 && j <= id_xsize)
			{
				int num_n = num_neighbours(grid, i, j);
				if (grid[i][j])
				{
					if (num_n != 2 && num_n != 3)
						new_grid[i][j] = false;
				}
				else if (num_n == 3) new_grid[i][j] = true;
			}
		}
	}
}


int main(int argc, char *argv[])
{

#ifdef TIME
	double tdif = -wTime();
#endif // TIME

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	srand(time(NULL) + id * 10);

	find_dimensions(p, proc_rows, proc_cols);	// Work out the best arrangement given the number of processors
	id_to_index(id, id_row, id_col);			// Work out the position of the proc in the grid of processors

	// Create main folder:
	string main_dir = to_string(proc_rows) + "x" + to_string(proc_cols) + "_" + to_string(max_steps);
	if (!filesystem::exists(main_dir))
	{
		filesystem::create_directory(main_dir);
	}

	// Work out the neighbors of the proc:
	vector<int> neighbors;						// [TL TM TR; ML __ MR; BL BM BR]
	int ix, jx;									
	for (int i = -1; i <= 1; i++)				// Start top left corner, end bottom right corner
	{
		for (int j = -1; j <= 1; j++)
		{
			if (i != 0 || j != 0)
			{
				if (periodic == 1)
				{
					ix = (i + id_row + proc_rows) % proc_rows;
					jx = (j + id_col + proc_cols) % proc_cols;

					neighbors.push_back(id_from_index(ix, jx));
				}
				else
				{
					ix = (i + id_row);
					jx = (j + id_col);
					
					if (ix > proc_rows || ix < 0 || jx > proc_cols || jx < 0)
					{
						//Out of bounds!
						neighbors.push_back(id);
					}
					else
					{
						neighbors.push_back(id_from_index(ix, jx));	
					}
				}
			}
		}
	}

	// Work out how many cells this proc should have:
	num_cells(jmax, id_col, id_xsize, proc_cols);
	num_cells(imax, id_row, id_ysize, proc_rows);

	//-----------------------
	// Prepare GAME OF LIFE:
	//-----------------------

	// Create new grids:
	bool **grid = new bool*[id_ysize + 2];
	bool **new_grid = new bool*[id_ysize + 2];

	// Initialise grids + paddings:
	for (int i = 0; i < id_ysize + 2; i++)
	{
		grid[i] = new bool[id_xsize + 2];
		new_grid[i] = new bool[id_xsize + 2];

		for (int j = 0; j < id_xsize + 2; j++)
		{
			new_grid[i][j] = 0;
			if (i >= 1 && i <= id_ysize && j >= 1 && j <= id_xsize)
			{
				grid[i][j] = (rand() % 2);
			}
			else grid[i][j] = 0;
		}
	}

	// Perform GAME OF LIFE:
	for (int n = 0; n < max_steps; n++)
	{
		MPI_Request* request1 = new MPI_Request[2 * 2];
		MPI_Request* request2 = new MPI_Request[6 * 2];

		bool *send_ML = new bool[id_ysize];
		bool *recv_ML = new bool[id_ysize];
		bool *send_MR = new bool[id_ysize];
		bool *recv_MR = new bool[id_ysize];

		for (int i = 0; i < id_ysize; i++)
		{
			send_ML[i] = grid[i + 1][1];
			send_MR[i] = grid[i + 1][id_xsize];
		}

		// SPECIFICALLY ADD ROW PADDING FIRST
		// TM -> BM	
		MPI_Isend(grid[1], id_xsize + 2, MPI_BYTE, neighbors[1], 2, MPI_COMM_WORLD, &request1[0]);
		MPI_Irecv(grid[id_ysize + 1], id_xsize + 2, MPI_BYTE, neighbors[6], 2, MPI_COMM_WORLD, &request1[1]);
		
		// BM -> TM
		MPI_Isend(grid[id_ysize], id_xsize + 2, MPI_BYTE, neighbors[6], 7, MPI_COMM_WORLD, &request1[2]);
		MPI_Irecv(grid[0], id_xsize + 2, MPI_BYTE, neighbors[1], 7, MPI_COMM_WORLD, &request1[3]);
		
		MPI_Waitall(4, request1, MPI_STATUSES_IGNORE);
		// THEN REPLACE THE PADDED CORNERS OF THE ROWS
		// TL -> BR
		MPI_Isend(&grid[1][1], 1, MPI_BYTE, neighbors[0], 1, MPI_COMM_WORLD, &request2[0]);
		MPI_Irecv(&grid[id_ysize + 1][id_xsize + 1], 1, MPI_BYTE, neighbors[7], 1, MPI_COMM_WORLD, &request2[1]);
		
		// TR -> BL
		MPI_Isend(&grid[1][id_xsize], 1, MPI_BYTE, neighbors[2], 3, MPI_COMM_WORLD, &request2[2]);
		MPI_Irecv(&grid[id_ysize + 1][0], 1, MPI_BYTE, neighbors[5], 3, MPI_COMM_WORLD, &request2[3]);
		
		// BL -> TR
		MPI_Isend(&grid[id_ysize][1], 1, MPI_BYTE, neighbors[5], 6, MPI_COMM_WORLD, &request2[4]);
		MPI_Irecv(&grid[0][id_xsize + 1], 1, MPI_BYTE, neighbors[2], 6, MPI_COMM_WORLD, &request2[5]);
		
		// BR -> TL
		MPI_Isend(&grid[id_ysize][id_xsize], 1, MPI_BYTE, neighbors[7], 8, MPI_COMM_WORLD, &request2[6]);
		MPI_Irecv(&grid[0][0], 1, MPI_BYTE, neighbors[0], 8, MPI_COMM_WORLD, &request2[7]);
		
		// ML -> MR
		MPI_Isend(send_ML, id_ysize, MPI_BYTE, neighbors[3], 4, MPI_COMM_WORLD, &request2[8]);
		MPI_Irecv(recv_ML, id_ysize, MPI_BYTE, neighbors[4], 4, MPI_COMM_WORLD, &request2[9]);
		
		// MR -> ML
		MPI_Isend(send_MR, id_ysize, MPI_BYTE, neighbors[4], 5, MPI_COMM_WORLD, &request2[10]);
		MPI_Irecv(recv_MR, id_ysize, MPI_BYTE, neighbors[3], 5, MPI_COMM_WORLD, &request2[11]);
		
		// WAIT FOR ALL COMMUNICATIONS
		MPI_Waitall(12, request2, MPI_STATUSES_IGNORE);

		for (int i = 0; i < id_ysize; i++)
		{
			grid[i + 1][id_xsize + 1] = recv_ML[i];
			grid[i + 1][0] = recv_MR[i];
		}

		delete[] send_ML;
		delete[] send_MR;
		delete[] recv_ML;
		delete[] recv_MR;

		do_iteration(grid, new_grid, id_xsize, id_ysize);
		swap(grid, new_grid);

		grid_to_file(n, main_dir, grid, id_xsize, id_ysize);
	}

#ifdef TIME
	tdif += wTime();									// Work out time difference.
	cout << "ID: " << id << "\tTime taken: " << tdif << endl;
#endif // TIME

	for (int i = 0; i < id_ysize; i++)
	{
		delete[] grid[i];
		delete[] new_grid[i];
	}
	delete[] grid;
	delete[] new_grid;

	MPI_Finalize();
	return 0;
}


//// check!
//stringstream stream2;
//stream2 << "ID: " << id << "\tTop+Bot+Corners+Left+Right" << endl;
//for (int i = 0; i < id_ysize + 2; i++)
//{
//	for (int j = 0; j < id_xsize + 2; j++)
//	{
//		stream2 << " " << grid[i][j];
//	}
//	stream2 << endl;
//}
//cout << stream2.str() << endl;