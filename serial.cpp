#include <iostream>
#include "math.h"
#include <fstream>
#include <ctime>
using namespace std;

int rows = 50, columns = 50;
int num_step = 100;

void periodic_neighbour(int *data)
{
    for (int j = 1; j < columns+1; j++)
    {
        data[0 * (columns+2) +j] = data[rows * (columns+2) +j];
        data[(rows+1) * (columns+2) +j] = data[1 * (columns+2) +j];
    }
    for (int i = 1; i < rows+1; i++)
    {
        data[i * (columns+2) + 0] = data[i * (columns+2) + columns];
        data[i * (columns+2) + (columns+1)] = data[i * (columns+2) + 1];
    }
    data[0] = data[rows * (columns+2) + columns];
    data[columns+1] = data[rows * (columns+2) + 1];
    data[(rows+1) * (columns+2)] = data[1 * (columns+2) + columns];
    data[(rows+1) * (columns+2) + (columns+1)] = data[1 * (columns+2) + 1];
}

void generate_randomly(int *data, bool periodic)
{
    for (int i = 0; i < rows+2; i++)
    {
        for (int j = 0; j < columns+2; j++)
            data[i * (columns+2) + j] = 0;
    }
    
    for (int i = 1; i < rows+1; i++)
    {
        for (int j = 1; j < columns+1; j++)
            data[i * (columns+2) + j] = rand()%2;
    }

    if (periodic == true)
        periodic_neighbour(data);
}

void generate_by_input(int *data, bool periodic)
{   
    for (int i = 0; i < rows+2; i++)
    {
        for (int j = 0; j < columns+2; j++)
            data[i * (columns+2) + j] = 0;
    }

    fstream fin;
    fin.open("3000.txt", fstream::in);
    if(!fin) {
        cout << "failed to open test file" << endl;
        return;
    }
    while (fin.good())
    {
        for (int i = 1; i < rows+1; i++)
        {
            for (int j = 1; j < columns+1; j++)
                fin >> data[i * (columns+2) + j];
        }
    }

    if (periodic == true)
    periodic_neighbour(data);
}


void print_data(int *data)
{
    for (int i = 1; i < rows+1; i++)
    {
        for (int j = 1; j < columns+1; j++)
            cout << data[i * (columns+2) + j] << "\t";
        cout << endl;
    }   
}

void print_data_with_layer(int *data)
{
    for (int i = 0; i < rows+2; i++)
    {
        for (int j = 0; j < columns+2; j++)
            cout << data[i * (columns+2) + j] << "\t";
        cout << endl;
    }   
}

void live_or_die(int *data1, int *data2, bool periodic)
{
    int time_step = 0;
    
    while (time_step < num_step)
    {
        for (int i = 1; i < rows + 1; i++)
        {
            for (int j = 1; j < columns + 1; j++)
            {
                int neighbour = 0;
                for (int k = i - 1; k <= i + 1; k++)
                {
                    for (int l = j - 1; l <= j + 1; l++)
                    {
                        neighbour += data1[k * (columns+2) + l];
                    }
                }
                neighbour -= data1[i * (columns+2) + j];
                if (data1[i * (columns+2) + j] == 1)
                {
                    if (neighbour == 2 || neighbour == 3) data2[i * (columns+2) + j] = 1;
                    else data2[i * (columns+2) + j] = 0;
                }
                else
                {
                    if (neighbour == 3) data2[i * (columns+2) + j] = 1;
                    else data2[i * (columns+2) + j] = 0;
                }                
            }
        }
        
        if (periodic)
            periodic_neighbour(data2);
        
        for (int i = 0; i < rows+2; i++)
        {
            for (int j = 0; j < columns+2; j++)
                data1[i * (columns+2) + j] = data2[i * (columns+2) + j];
        }
        // data_old = data_new;
        time_step++;
    }
}


int main()
{   
    int *data_old = new int[(rows+2) * (columns+2)];
    int *data_new = new int[(rows+2) * (columns+2)];
    // generate_by_input(data_old, true);
    generate_randomly(data_old, true);

    // cout << "before:" << endl;
    // print_data(data_old);

    clock_t startTime, endTime;
    startTime = clock();

    live_or_die(data_old, data_new, true);

    endTime = clock();
    double running_time = (double) (endTime - startTime) / CLOCKS_PER_SEC;
    cout << "The running time is: " << running_time << "s" << endl;

    // cout << "after:" << endl;
    // print_data(data_old);
    return 0;
}