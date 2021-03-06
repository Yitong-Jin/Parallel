# ACSE-6 individual assignment - Game of Life

* Name: Yitong Jin

* CID: 01776029

* Email: yj319@ic.ac.uk

## Overview

This work is aim to write a parallel version of Conway’s Game of Life by using c-plus-plus and MPI. This assignment includes several parts which are as follows:

* readme file of this assignment (the file you are reading now): "README.md"

* serial version code of Game of Conway’s Life: "serial.cpp"

* parallel version code of Game of Conway’s Life: "parallel.cpp"

* directory for storing output data: "output_data/"

* directory for test examples and solutions: "test_example/"

* three txt files that I used for performance evaluation: "500.txt", "3000.txt" and "5000.txt" 

* a xlsx file generated by EXCEL which contains the data and figures of performance: "performance_data.xlsx"

* postprocessing script to merge separate generated by each processor and to visualize simulation of the game: "Postprocessing.ipynb"

* a script of makefile for HPC: "Makefile"

* a script of setting environment of HPC: "my_script.pbs"

* a report for analysis and discussion of the program’s performance: "performance_report.pdf"

## Environment

* operating system: MacOS Mojave 10.14.5

* Programming language: c++

* Compiler: mpic++

## Compile

For parallel version code:
> mpic++ parallel.cpp -o parallel

For serial version code:
> g++ serial.cpp -o serial

## Run

Parallel code: (here we use 4 processor to run the program as example)
> mpirun -np 4 ./parallel

Serial code:
> ./serial

## User guidance

In this assignment the serial version code is just a reference for checking speedup ratio by running parallel version code. You can run serial code yourself or just ignore it.

For parallel version code, there are some essential parameters which determine the game status, namely:
* **bool** periodic; (true for periodic boundary, false for fixed boundary)
* **int** rows; (the number of rows in the domain)
* **int** columns; (the number of columns in the domain)
* **int** num_step; (the total number of steps of the game)
* **bool** input_randomly; (true for generate intial configuration, false for read data from txt file)
* **char** input_file[30]; (only useful when "input_randomly" = false, enter the file path of the txt file which contains the inital data)
* **bool** data_output; (true if you want data output, otherwise choose false)

I have set default values for above parameters, which will result in a Conway’s Game of Life with a 50 x 50 domain size, the total number of step is 100, input initial condition randomly, output data will be created in the directory named "output_data". 

### **IMPORTANT**
If you want to running the program with different game setting, you should change above 7 parameters before running the program. And if you set variable data_output as true that means you will get the output data of every single processor at every step. To avoid confusion, you need to delete all of data files already in "output_data/" directory before you run the program.

## Introduction of code and process

* START
* Before running program, setting parameters and empty "output_data" directory.
* Compile and run "parallel.cpp".
* Initial configuration generated by random or read from txt file.
* The program decompose the whole domain by several horizontal stripes, as many as the total number of processors we use.
* According to total number of processors we use, allocate appropriate number of rows to each processor. Also allocate two more rows to each processor as the up padding and down padding in order to receive the data from adjacent processors. There is no need to append left and right paddings because the neighbours' information of the first column and the last column can be obtained from same processor. I just use a specific algorithm to caculate their alive neighbours. Named this space as "data_old".
* Creat one more same size memory for each processor to store the update living conditions for each cell in every single update step. Named this space as "data_new".
* Assign intial data to each processor storing in the memory whcih just allocate to each processor with the up padding and down padding blank.
* Use MPI to let each processor communicate with their adjacent processors to send the first row and the last row to the previous and next processors respectively. Meanwhile, receive the last row of previous processor storing in the up padding and the first row of next processor storing in the down padding.
* Now the Game of Life starts and timing starts here.
* Update for "num_step" times to get the final result. In each processor, loop from the first cell to the last cell to caculate the number of alive neighbors for each cell. Find out the next step living condition for every cell and update living condition by storing in "data_new".
* After every single time update of living condition, each processor communicates with their adjacent two processors to swap the information of thier first and last row, just as same as what they do in first time communication.
* Assign the data in "data_new" to "data_old" before next update of living condition. As the "data_new" in previous step is exactly same as the "data_old" of next step.
* Cumulate the number of loop until the step equal to total steps.
* Timing end and compare the running time between every processors, pick the longest running time as the program running time.
* Close "parallel.cpp" and you can find the output data in "output_data/" directory.
* Open "Postprocessing.ipynb" by jupyter notebook, run the first code cell, merged data will also generate in "output_data/" directory. Merged data concatenates data from every processors in every single step.
* Run the second code cell, simulation animation can be generated.
* END

*More detailed comments and introduction can be find in both "parallel.cpp" and "Postprocessing.ipynb".*







