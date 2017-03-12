#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include <time.h>

#define COLS 10000
#define ROWS 10000
#define TEMP 50.
#define DEBUG 0
#define EPS 1e-3
#define I_FIX 2
#define J_FIX 2

double max_abs(double** m1, double** m2, int row){
    double max_val = DBL_MIN;
    for (int i = 1; i < row+1; i++)
        for (int j = 0; j < COLS+2; j++){
            if (fabs(m1[i][j] - m2[i][j]) > max_val) {
                max_val = fabs(m1[i][j] - m2[i][j]);
            }
        }
    return max_val;
}

double get_max(double* m, int p)
{
	double max = DBL_MIN;
	for(int i =0; i<p; i++)	
		if(m[i] > max)
			max = m[i];

	return max; 
}

void copy_top_ghost(double *temp, double **dest_matrix)
{
	for (int i=1; i < COLS+1; i++)
		dest_matrix[0][i] = temp[i];
}


void copy_bottom_ghost(double *temp, double **dest_matrix, int w)
{
	for (int i=1; i < COLS+1; i++)
		dest_matrix[w+1][i] = temp[i];
}

void print_matrix(double** matrix, int w){
    for (int i = 1; i < w+1; i++) {
        for (int j = 1; j < COLS+1; j++)
            printf("%f ", matrix[i][j]);
        printf("\n");
    }
}

void copy_matrix(double** dest, double** source, int w) {
    for (int i = 0; i < w+2; i++)
        for (int j = 0; j < COLS+2; j++)
            dest[i][j] = source[i][j];
}

double** alloc_matrix(int w){
    double** matrix;
    matrix = (double**) malloc((w+2) * sizeof(double *));
    matrix[0] = (double*) malloc((w+2) * (COLS+2) * sizeof(double));
    for (int i = 1; i < (w+2); i++)
        matrix[i] = matrix[0] + i*(COLS+2);
    return matrix;
}

void compute_new_values(double** old_matrix, double** new_matrix, int row, int my_rank, int p ){
    int h = ROWS/p;
    for (int i = 1; i < row+1; i++)
        for (int j= 1; j < COLS+1; j++)
            new_matrix[i][j] =
                    0.25 * (old_matrix[i-1][j] + old_matrix[i+1][j] +
                            old_matrix[i][j-1] + old_matrix[i][j+1]);

    if ((my_rank * h)  <= I_FIX && I_FIX < (my_rank+1)*row )
    	new_matrix[I_FIX - (my_rank*h) +1][J_FIX] = TEMP;

    //new_matrix[I_FIX][J_FIX] = TEMP;
}

void init_matrix(double** matrix, int w, int my_rank, int p){
    for (int i = 0; i < w+2; i++)
        for (int j = 0; j < COLS+2; j++) {
            matrix[i][j] = 0.;
        }
    int h = (int) (ROWS/p);
    if ((my_rank * h)  <= I_FIX && I_FIX < (my_rank+1)*w )
    	matrix[I_FIX - (my_rank*h)+1 ][J_FIX] = TEMP;
}

void update_top_ghost(double** matrix)
{
	for( int j=1; j < COLS+1; j++)
		matrix[0][j] = matrix[1][j];
}

void update_bottom_ghost(double** matrix, int w)
{
	for( int j=1; j < COLS+1; j++)
		matrix[w+1][j] = matrix[w][j];
}

int get_num_rows(int my_rank, int p)
{
	int r, h;
	h = (int) (ROWS/p); 
	r = ROWS%p;

	if(my_rank == p-1)
		return h+r;
	else
		return h;
}

int main(int argc, char *argv[]) {

	int p, my_rank, my_top, my_bottom;
	int w; //#of rows per process
	int tag=0;
	double l_max_diff;
	double start, end; //To measure run time

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&p);
	start = MPI_Wtime();
	if ( p > ROWS)
	{
		printf("Too many number of processes");
		exit(0);
	}

	if ( I_FIX < 0 || I_FIX > ROWS)
	{
		printf("Invalid arguments for heat source");
		exit(0);
	}
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

	
	MPI_Status status;
	/*initialize local matrix and parameters*/
	w = get_num_rows(my_rank, p);
	if ( w < 0 )
	{
		printf("Error");
		exit(0);
	}
	my_top = my_rank-1;
	my_bottom = my_rank+1;


    	double **a_old = alloc_matrix(w); //allocate memory for local the matrices
    	double **a_new = alloc_matrix(w);

    	init_matrix(a_old, w, my_rank, p); //initialize the matrices
    	init_matrix(a_new, w, my_rank, p);

	double **final;
//	if (my_rank == 0)
//	{
		 final = alloc_matrix(ROWS);
//	}
	if (DEBUG) {
	    printf("w = %d, my_top = %d, my_bottom = %d, my_rank = %d, #pr = %d \n", w, my_top, my_bottom, my_rank, p);
        }

	
	MPI_Request top_recv, bottom_send;
	MPI_Request top_send, bottom_recv;
	double *temp_bottom;
	double *temp_top;
	if (my_rank > 0)
	{
		temp_top = malloc(sizeof(double) * (COLS+2));
	}
	if(my_rank < p-1)
	{
		temp_bottom = malloc(sizeof(double) * (COLS+2));
	}
	
	
    while (1) {
     	//MPI_Barrier(MPI_COMM_WORLD); 
	
        if (DEBUG)
            printf("Performing a new iteration...\n");

        //compute new values and put them into a_new

        compute_new_values(a_old, a_new, w, my_rank, p);

        if (DEBUG) {
           // printf("a_old at p = %d is:\n",my_rank); //output matrix to screen
           // print_matrix(a_old, w);

           // printf("a_new at p = %d is:\n",my_rank);
           // print_matrix(a_new, w);
        }

        //calculate the maximum absolute differences among pairwise
        // differences of old and new matrix elements
            l_max_diff = max_abs(a_old, a_new, w);

	double *all_max_diff = malloc(sizeof(double)*p);

	MPI_Allgather(&l_max_diff, 1, MPI_DOUBLE, all_max_diff, 1, MPI_DOUBLE, MPI_COMM_WORLD);

	double g_max_diff = get_max(all_max_diff, p);
	
        if (DEBUG)
	{
            printf("Max local diff at p = %d is: %f\n",my_rank, l_max_diff);
	    printf("Max global diff at p = %d is: %f\n", my_rank, g_max_diff);
	}

	if (g_max_diff < EPS)
		break;

	copy_matrix(a_old, a_new, w); //assign values of a_new to a_old
	/* send top and receive bottom ghost cell */

	if (my_rank > 0)
	{
		update_top_ghost(a_new);
		MPI_Isend(a_new[0], COLS+2, MPI_DOUBLE, my_top, tag, MPI_COMM_WORLD, &top_send);
		MPI_Irecv(temp_top, COLS+2,MPI_DOUBLE, my_top, tag, MPI_COMM_WORLD, &bottom_recv);

		
	}

	if (my_rank < p-1)
	{
		update_bottom_ghost(a_new, w);
		MPI_Isend(a_new[w+1], COLS+2, MPI_DOUBLE, my_bottom, tag, MPI_COMM_WORLD, &bottom_send);
		MPI_Irecv(temp_bottom, COLS+2,MPI_DOUBLE, my_bottom, tag, MPI_COMM_WORLD, &top_recv);

	}
	//top
	if (my_rank > 0)
	{
		MPI_Wait(&bottom_recv,&status);
		//copy_top_ghost(temp_top, a_old);
		copy_top_ghost(temp_top, a_new);
	}
		
	//bottom
	if (my_rank < p-1)
	{
		MPI_Wait(&top_recv,&status);
		//copy_bottom_ghost(temp_bottom, a_old, w);
		copy_bottom_ghost(temp_bottom, a_new, w);
	}
	
	if (my_rank > 0)
	{
		MPI_Wait(&top_send,&status);
	}
	if (my_rank < p-1)
	{
		MPI_Wait(&bottom_send,&status);
	}
	//MPI_Barrier( MPI_COMM_WORLD );
	copy_matrix(a_old, a_new, w); //assign values of a_new to a_old
        if (DEBUG)
            printf("End of iteration\n\n");
    }
	if (DEBUG)
	{
		printf("\nExited Loop, p = %d\n",my_rank);
	}
 
//send metrics to root

	MPI_Gather( a_new[1], (ROWS/p)*(COLS+2), MPI_DOUBLE, final[1],(ROWS/p)*(COLS+2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (DEBUG)
	{
		printf("\nGather Successful, p=%d\n",my_rank);
	}
	if (my_rank == p-1)
	{

		if (DEBUG)
		{
			printf("\n sending extra\n");
		}
		if ( p * w > ROWS)
		{
			int offset = 1+ (int)(ROWS/p);
			MPI_Send(a_new[offset],(w-offset+1)*(COLS+2), MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
		}
	
	}

	if (my_rank == 0)
	{
		if (DEBUG)
		{
			printf("\n receiving extra\n");
		}
		if ( p*w < ROWS)
			MPI_Recv(final[w*p+1],(ROWS-(w*p))*(COLS+2), MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
		if(DEBUG)
		{
			printf("Root received all");
		}
		end = MPI_Wtime();
    		print_matrix(final,ROWS);
		printf("\nTotal Execution Time = %f\n", end-start);
	}
	MPI_Finalize();
    return 0;
}
