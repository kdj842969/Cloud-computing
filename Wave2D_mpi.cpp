#include <iostream>
#include "Timer.h"
#include <stdlib.h>   // atoi
#include "mpi.h"
#include <omp.h>

int default_size = 100;  // the default system size
int defaultCellWidth = 8;
double c = 1.0;      // wave speed
double dt = 0.1;     // time quantum
double dd = 2.0;     // change in system
double cons1 = c*c*(dt/dd)*(dt/dd); //constant value 1
double cons2 = c*c*(dt/dd)*(dt/dd)/2; //constant value 2

using namespace std;

int main( int argc, char *argv[] ) {
    // verify arguments
    if ( argc != 5 ) {
        cerr << "usage: Wave2D size max_time interval" << endl;
        return -1;
    }
    int size = atoi( argv[1] );
    int max_time = atoi( argv[2] );
    int interval  = atoi( argv[3] );
    int my_rank = 0;            // used by MPI
    int mpi_size = 1;           // used by MPI
    MPI_Status status;
    int threads = atoi(argv[4]);
    
    if ( size < 100 || max_time < 3 || interval < 0 ) {
        cerr << "usage: Wave2D size max_time interval" << endl;
        cerr << "       where size >= 100 && time >= 3 && interval >= 0" << endl;
        return -1;
    }
    
    MPI_Init(&argc, &argv); // start MPI
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );

    omp_set_num_threads(threads);
    // create a simulation space
    double z[3][size][size];
    for ( int p = 0; p < 3; p++ )
        for ( int i = 0; i < size; i++ )
            for ( int j = 0; j < size; j++ )
                z[p][i][j] = 0.0; // no wave
    
    // start a timer
    Timer time;
    time.start( );
    
    // time = 0;
    // initialize the simulation space: calculate z[0][][]
    int weight = size / default_size;
#pragma omp parallel for shared(z)
    for( int i = 0; i < size; i++ ) {
        for( int j = 0; j < size; j++ ) {
            if( i > 40 * weight && i < 60 * weight  &&
               j > 40 * weight && j < 60 * weight ) {
                z[0][i][j] = 20.0;
            } else {
                z[0][i][j] = 0.0;
            }
        }
    }
    
#pragma omp parallel shared(z)
    for(int i = 0; i < size; i++ ) {
        for( int j = 0; j < size; j++ ) {
            if((i==0)||(i==(size-1))||(j==0)||(j==(size-1))) {
                z[1][i][j]=0.0;
            } else {
                z[1][i][j]=z[0][i][j]+(c*c/2)*(dt/dd)*(dt/dd)*(z[0][i+1][j]+z[0][i-1][j]+z[0][i][j+1]+z[0][i][j-1]-4.0*z[0][i][j]);
            }
        }
    }

    // time = 1
    // calculate z[1][][]
    // cells not on edge
    // IMPLEMENT BY YOURSELF !!!
    
    int stripe = size / mpi_size;     // partitioned stripe stripe this array into 4 stripes, mpi_size = 4
    // if size=400, mpi_size=4, stripe_size = 100; pass 400, 800, no remainer
    
    
    //if(my_rank==0){
        // master sends each partition of a[] to a different slave
        // master also sends b[] to all slaves
      //  for(int i=1; i<mpi_size;i++){
        //    MPI_Send(z, 3*size*size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD); //a must point to different stripe
        //}
    //} else {
        //slaves receive
      //  MPI_Recv(z, 3*size*size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    //}
    
    
    // simulate wave diffusion from time = 2
    for ( int t = 2; t < max_time; t++ ) {
#pragma omp parallel shared(z)
        for(int i = my_rank*stripe; i<(my_rank+1)*stripe; i++){
            for(int j=0; j < size; j++){
                if((i==0)||(i==(size-1))||(j==0)||(j==(size-1))) {
                    z[t%3][i][j]=0.0;
                } else {
                    z[t%3][i][j]=2.0*z[(t-1)%3][i][j]-z[(t-2)%3][i][j]+c*c*(dt/dd)*(dt/dd)*(z[(t-1)%3][i+1][j]+z[(t-1)%3][i-1][j]+z[(t-1)%3][i][j+1]+z[(t-1)%3][i][j-1]-4.0*z[(t-1)%3][i][j]);
                }
            }
        }

        int left = my_rank-1;
        int right = my_rank+1;
        if (left < 0) {
            left = MPI_PROC_NULL;
        }
        if (right == mpi_size) {
            right = MPI_PROC_NULL;
        }
         //cout<<"rank="<<my_rank<<endl;
        if(my_rank%2 == 0){
            MPI_Send(&z[t%3][my_rank*stripe], size, MPI_DOUBLE, left, 0, MPI_COMM_WORLD);
            MPI_Send(&z[t%3][(my_rank+1)*stripe-1], size, MPI_DOUBLE, right, 0, MPI_COMM_WORLD);
            MPI_Recv(&z[t%3][my_rank*stripe-1], size, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&z[t%3][(my_rank+1)*stripe], size, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
        } else {
            MPI_Recv(&z[t%3][my_rank*stripe-1], size, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&z[t%3][(my_rank+1)*stripe], size, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &status);
            MPI_Send(&z[t%3][my_rank*stripe], size, MPI_DOUBLE, left, 0, MPI_COMM_WORLD);
            MPI_Send(&z[t%3][(my_rank+1)*stripe-1], size, MPI_DOUBLE, right, 0, MPI_COMM_WORLD);
        }
        
	if((t%interval)==0){
            if(my_rank == 0){
                for(int i = 1; i<mpi_size;i++){
                    MPI_Recv(&z[t%3][i*stripe], stripe*size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status); //receive results, to appropriate partition
                }
            } else {
                MPI_Send(&z[t%3][my_rank*stripe], stripe*size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }  
            if(my_rank==0){
                cout<<t<<endl;
                for(int i=0; i<size;i++){
                    for(int j=0; j<size;j++){
                        cout<<z[t%3][i][j]<<" ";
                    }
                    cout<<endl;
                }
            }
        }
        
    } // end of simulation
    

    if(my_rank == 0){
        for(int i = 1; i<mpi_size;i++){
            MPI_Recv(&z[(max_time-1)%3][i*stripe], stripe*size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status); //receive results, to appropriate partition
        }
    } else {
           MPI_Send(&z[(max_time-1)%3][my_rank*stripe], stripe*size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
  
                         
    if(my_rank==0){
        cout<<max_time-1<<endl;
        for(int i=0; i<size;i++){
            for(int j=0; j<size;j++){
                cout<<z[(max_time-1)%3][i][j]<<" ";
            }
            cout<<endl;
        }
    }
    
    


    // finish the timer
    cerr << "Elapsed time = " << time.lap( ) << endl;
    MPI_Finalize();
}

