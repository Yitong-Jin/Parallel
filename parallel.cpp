#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
using namespace std;

// Some functions in this program, definiton and commenting after main function
void generate_by_input(int *data);
void generate_randomly(int *data);
void update(int *data1, int *data2, int num_rows, int num_columns);
void allocate_num_rows(int *start_row, int *end_row);
void allocate_data(int num_rows, int num_columns, int *data_all, int *data_old);
void communication(int *data1, int num_rows, int num_columns, bool periodic);
void output_file(int *data, int num_rows, int num_columns, int time_step);

// Some global variables
int id, p;
int tag_num = 1;
int *data_all = nullptr;
int *start_row, *end_row;
int *data_old = nullptr;
int *data_new = nullptr;
int *data_to_send = nullptr;
int *data_to_recv = nullptr;

// IMPORTANT ! READ BEFORE RUNNING THE PROGRAM (also mentioned in README file)
// Here are some initial conditions as following, you may need to 
// set them yourself according to specific game initial condition.

// you need to switch "periodic" from false to true if you want a periodic game domain.
bool periodic = true;

// set game domain before the game
int rows = 50, columns = 50;

// set total steps of game before running the progarm
int num_step = 100;

// if you want to generate initial configuration randomly, switch "input_randomly" to "true".
// OR if you want to input initial configuration from a "txt" file, switch "input_randomly"
// to "false" and assign the filepath to variable "input_file".
bool input_randomly = true;
char input_file[30] = "500.txt";

// if you want to output data file, switch "data_output" to file and you will find output files
// in directory "./output_data/". Again and again, please clear the "./output_data/" directory 
// before you want to have new output data files. Also add a blank "./output_data/" directory when
// you are running HPC if you want to output data after running by HPC. For merge data and visualization
// please run the code cell in "Postprocessing.ipynb".
bool data_output = true;

int main(int argc, char *argv[])
{   
    // generate initial distributions
    data_all = new int[rows * columns];
    if(input_randomly == true) generate_randomly(data_all);
    else generate_by_input(data_all);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	srand(time(NULL) + id * 1000);

    // allocate the number of row to each processor
    start_row = new int[p];
    end_row = new int[p];
    allocate_num_rows(start_row, end_row);

    // for each processor, allocates enough memory to store its portion of the domain
    // and the same size memory to store the next step data to finish one step update
    // here each processor is assigned two more rows to receive the data from the last
    // row of previous processor and the first row of next processor
    data_old = new int[(end_row[id] - start_row[id] + 1 + 2) * columns];
    data_new = new int[(end_row[id] - start_row[id] + 1 + 2) * columns];
    int num_rows = end_row[id] - start_row[id] + 1 + 2;
    int num_columns = columns;
    // allocate initial data to each processor
    allocate_data(num_rows, num_columns, data_all, data_old);
    
    // at initial state, send the first row to previous processor and the last row to next processor.
    // Also receive the first row of next processor and the last row of previous processor to get
    // necessary information to finish every step update.
    communication(data_old, num_rows, num_columns, periodic);
    // output initial status
    if (data_output == true)
        output_file(data_old, num_rows, num_columns, 0);

    // timing start
    clock_t startTime, endTime;
    startTime = clock();

    // the loop for update
    int time_step = 0;
    while (time_step < num_step)
    {
        // caculating number of alive neighbour and finding out survival status for every grid in domain
        update(data_old, data_new, num_rows, num_columns);
        time_step++;
        //output the data for every step of every processor
        if (data_output == true)
            output_file(data_new, num_rows, num_columns, time_step);
    }
    
    // timing end
    endTime = clock();
    double running_time_each = (double) (endTime - startTime) / CLOCKS_PER_SEC;
    // cout << "The running time of processor " << id << " is: " << running_time_each << "s" << endl;
    double running_time;
    // each processor may have different running time, here we need to take the longest time among 
    // all of processors, because only every processors finish their jobs can we get the final result. 
    MPI_Reduce(&running_time_each, &running_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (id == 0) cout << "The running time is: " << running_time << " s." << endl;

    // if we called new, we call delete!
    delete[] data_all;
    delete[] start_row;
    delete[] end_row;
    delete[] data_old;
    delete[] data_new;
    cout.flush();
    MPI_Finalize();
    return 0;
}



/** Generate initial configuration by input
        
    read the data in specific txt file one by one and allocate to array

    @param data       the pointer of data array to receive data from the txt file
*/
void generate_by_input(int *data)
{   
    fstream fin;
    fin.open(input_file, fstream::in);
    if(!fin) {
        cout << "failed to open test file" << endl;
        return;
    }
    while (fin.good())
    {        
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
                fin >> data[i * columns + j];
        }      
    }
    fin.close();
}


/** Generate initial configuration by random
        
    fullfill the array by random number to construct initial state

    @param data       the pointer of data array to receive random number
*/
void generate_randomly(int *data)
{   
    srand(time(NULL) + id * 1000);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
            data[i * columns + j] = rand()%2;
    }
}


/** Allocate appropriate number of rows to each processor
        
    construct two list start_row and end_row and the size of each equal to the number of processor
    for the whole domain, we want to divide it into "p"(the number of processor) as evenly as possible
    by calling this function, we can store the starting row of i th processor at start_row[i] and store
    the ending row of i th processor at end_row[i]. For example, the whole domain has 8 rows, and we have 3
    processor available. Now we can divide 8 rows as 2, 3, 3 and assign them to processor 0, 1, 2 respectively.
    In this situation, start_row =[0,2,5] and end_row = [1,4,7]. 

    @param start_row     the pointer of array to store the serial number of starting row of each processor.
    @param end_row       the pointer of array to store the serial number of ending row of each processor.
*/
void allocate_num_rows(int *start_row, int *end_row)
{
    int row_rem = rows;
    start_row[0] = 0;

    for (int i = 0; i < p; i++)
    {
        int row_on_proc = row_rem / (p - i);
        end_row[i] = start_row[i] + row_on_proc - 1;

        if (i + 1 < p)
            start_row[i + 1] = end_row[i] + 1;
        row_rem -= row_on_proc;
    }
}


/** Allocate data to each processor

    @param num_rows         number of rows that assign to each processor
    @param num_columns      number of columns that assign to each processor
    @param data_all         the pointer of data for the whole game domain
    @param data_old         to receive the appropriate portion of data sliced from data_all, also the data before every single update
*/
void allocate_data(int num_rows, int num_columns, int *data_all, int *data_old)
{
    for (int i = 1; i < num_rows - 1; i++)
    {
        for (int j = 0; j < num_columns; j++)
            data_old[i * num_columns + j] = data_all[(start_row[id] + i - 1) * columns + j];
    }
}


/** Communicate to adjacent processor
        
    Processors send the first row to previous processor and the last row to next processor
    Meanwhile receive the first row of next processor and the last row of previous processor
    However, for non-periodic domain, the 0 processor do not have previous processor and the
    "p-1"th processor do not have next processor, so some communication should be omit.

    @param data1            the pointer of data for each processor
    @param num_rows         number of rows that assign to each processor
    @param num_columns      number of columns that assign to each processor
    @param periodic         game domain is periodic or not
*/
void communication(int *data1, int num_rows, int num_columns, bool periodic)
{
    MPI_Request request;
    if (id == 0 && periodic == false)
    {
        for (int j = 0; j < num_columns; j++)
            data1[j] = 0;
        MPI_Isend(&data1[(num_rows - 2) * num_columns], num_columns, MPI_INT, 1, tag_num, MPI_COMM_WORLD, &request);
        MPI_Recv(&data1[(num_rows - 1) * num_columns], num_columns, MPI_INT, 1, tag_num, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else if (id == p-1 && periodic == false)
    {
        for (int j = 0; j < num_columns; j++)
            data1[(num_rows-1) * num_columns + j] = 0;
        MPI_Isend(&data1[num_columns], num_columns, MPI_INT, p-2, tag_num, MPI_COMM_WORLD, &request);
        MPI_Recv(&data1[0], num_columns, MPI_INT, p-2, tag_num, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        MPI_Isend(&data1[(num_rows - 2) * num_columns], num_columns, MPI_INT, (id + 1) % p, tag_num, MPI_COMM_WORLD, &request);
        MPI_Isend(&data1[num_columns], num_columns, MPI_INT, (id + p - 1) % p, tag_num, MPI_COMM_WORLD, &request);        
        MPI_Recv(&data1[(num_rows - 1) * num_columns], num_columns, MPI_INT, (id + 1) % p, tag_num, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&data1[0], num_columns, MPI_INT, (id + p - 1) % p, tag_num, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}


/** Counting alive neighbour and find out next step situation
        
    caculating alive neighbour for every grid in domain, get next step living information
    communicating with adjacent processor to update the information of extra padding 
    assign new data for previous step to old data for next step

    @param data_old         the pointer of data befor every single update
    @param data_new         the pointer of data after every single update
    @param num_rows         number of rows that assign to each processor
    @param num_columns      number of columns that assign to each processor
*/
void update(int *data_old, int *data_new, int num_rows, int num_columns)
{
    // counting alive neighbour and find out live or dead at next step
    for (int i = 1; i < num_rows - 1; i++)
    {
        for (int j = 0; j < num_columns; j++)
        {
            int neighbour = 0;
            for (int k = i - 1; k <= i + 1; k++)
            {
                for (int l = j - 1; l <= j + 1; l++)
                {
                    if (periodic == false && (l < 0 || l > num_columns - 1)) continue;
                    else if (periodic == true && l < 0) 
                        neighbour += data_old[k * num_columns + num_columns - 1];
                    else if (periodic == true && l > num_columns - 1)
                        neighbour += data_old[k * num_columns];
                    else neighbour += data_old[k * num_columns + l];       
                }
            }
            neighbour -= data_old[i * num_columns + j];
            if (data_old[i * num_columns + j] == 1)
            {
                if (neighbour == 2 || neighbour == 3) data_new[i * num_columns + j] = 1;
                else data_new[i * num_columns + j] = 0;
            }
            else
            {
                if (neighbour == 3) data_new[i * num_columns + j] = 1;
                else data_new[i * num_columns + j] = 0;
            }               
        }
    }

    // swaping the information in up and down padding after next step neighbour counting
    communication(data_new, num_rows, num_columns, periodic);

    // before next update operation, set the new data for previous step as the old data for next step      
    for (int i = 0; i < num_rows; i++)
    {
        for (int j = 0; j < num_columns; j++)
            data_old[i * num_columns + j] = data_new[i * num_columns + j];
    }
}


/** Generate output file at each step each processor
        
    At each step each process, generate their own output file. The file name of output file includes
    the total processors we use, total steps, id of specific proessor and specific step. All of files
    generate will store in the directory "output_data".

    @param data             the pointer of data to be output
    @param num_rows         number of rows that assign to each processor
    @param num_columns      number of columns that assign to each processor
    @param time_step        the current time_step
*/
void output_file(int *data, int num_rows, int num_columns, int time_step)
{
    char output_file[50];
    sprintf(output_file, "output_data/p%ds%d_%d,%d.txt", p, num_step, id, time_step);
    fstream fout;
    fout.open(output_file, ios::out);
    if(!fout) {
        cout << "failed to write output file" << endl;
        return;
    }
    for (int i = 1; i < num_rows - 1; i++)
    {
        for (int j = 0; j < num_columns; j++)
            fout << data[i * num_columns + j];
        fout << endl;
    }
    fout.close();
}


