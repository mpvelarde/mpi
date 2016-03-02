#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "stack.h"

#define EPSILON 1e-3
#define F(arg)  cosh(arg)*cosh(arg)*cosh(arg)*cosh(arg)
#define A 0.0
#define B 5.0

#define SLEEPTIME 1

#define MAX_TASKS 100
#define NO_MORE_TASKS MAX_TASKS+1

int *tasks_per_process;

double farmer(int);
void worker(int);
double quad (double, double, double, double, double);

int main(int argc, char **argv ) {
  int i, myid, numprocs;
  double area, a, b;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if(numprocs < 2) {
    fprintf(stderr, "ERROR: Must have at least 2 processes to run\n");
    MPI_Finalize();
    exit(1);
  }

  if (myid == 0) { // Farmer
    // init counters
    tasks_per_process = (int *) malloc(sizeof(int)*(numprocs));
    for (i=0; i<numprocs; i++) {
      tasks_per_process[i]=0;
    }
  }

  if (myid == 0) { // Farmer
    area = farmer(numprocs);
  } else { //Workers
    worker(myid);
  }

  if(myid == 0) {
    fprintf(stdout, "Area=%lf\n", area);
    fprintf(stdout, "\nTasks Per Process\n");
    for (i=0; i<numprocs; i++) {
      fprintf(stdout, "%d\t", i);
    }
    fprintf(stdout, "\n");
    for (i=0; i<numprocs; i++) {
      fprintf(stdout, "%d\t", tasks_per_process[i]);
    }
    fprintf(stdout, "\n");
    free(tasks_per_process);
  }
  MPI_Finalize();
  return 0;
}


double farmer(int numprocs) {
    // Init stack of tasks
    stack *tasks = new_stack();
    
    // Init variables to receive info from MPI_Recv
    int temp, tag, who;
    MPI_Status status;
    double area[MAX_TASKS], result[MAX_TASKS];
    double points[5];
    double *task;
    
    // Generate tasks
    int i = 0;
    while(i < MAX_TASKS) {
        points[0] = A;
        points[1] = B;
        points[2] = F(A);
        points[3] = F(B);
        points[4] = (F(A)+F(B)) * (B-A)/2;
        
        push(points, tasks);
        i++;
    }
    
    // Assume at least as many tasks as workers
    for (i=0; i < (numprocs-1); i++) {
        task = pop(tasks);
        MPI_Send(&task, 5, MPI_DOUBLE, i+1, i, MPI_COMM_WORLD);
    }
    
    while (i<MAX_TASKS) {
        MPI_Recv(&temp, 2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        who = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        result[tag] = temp;
        
        task = pop(tasks);
        MPI_Send(&task, 5, MPI_DOUBLE, who, i, MPI_COMM_WORLD);
        i++;
    }
    
    for (i=0; i < numprocs; i++) {
        MPI_Recv(&temp, 2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        who = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        result[tag] = temp;
        
        task = pop(tasks);
        MPI_Send(&task, 1, MPI_DOUBLE, who, NO_MORE_TASKS, MPI_COMM_WORLD);
    }
    printf("Farmer for %d wuu \n", numprocs);
}

void worker(int mypid) {
    // Init counters for worker
    int tasksdone = 0;
    double workdone = 0;
    
    // Init variables to receive task
    double task[5];
    int tag;
    double result[2];
    double left, right, fleft, fright, lrarea;
    double mid, fmid, larea, rarea;
    MPI_Status status;
    
    // Receive task
    MPI_Recv(&task, 5, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    tag = status.MPI_TAG;
    
    while (tag != NO_MORE_TASKS) {
        // Get variables
        left = task[0];
        right = task[1];
        fleft = task[2];
        fright = task[3];
        lrarea = task[4];
        
        mid = (left + right) / 2;
        fmid = F(mid);
        
        // Get result
        result[0] = (fleft + fmid) * (mid - left) / 2;
        result[1] = (fmid + fright) * (right - mid) / 2;
        
        // Update counters
        workdone+= result[0] + result[1];
        tasksdone++;
        
        // Send result
        MPI_Send(&result, 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        
        // Receive next task
        MPI_Recv(&task, 5, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        tag = status.MPI_TAG;
    }
    printf("Worker %d solved %d tasks totalling %f units of work \n", mypid, tasksdone, workdone);
}

double quad(double left, double right, double fleft, double fright, double lrarea) {
    double mid, fmid, larea, rarea;
    
    mid = (left + right) / 2;
    fmid = F(mid);
    larea = (fleft + fmid) * (mid - left) / 2;
    rarea = (fmid + fright) * (right - mid) / 2;
    if( fabs((larea + rarea) - lrarea) > EPSILON ) {
        larea = quad(left, mid, fleft, fmid, larea);
        rarea = quad(mid, right, fmid, fright, rarea);
    }
    return (larea + rarea);
}


