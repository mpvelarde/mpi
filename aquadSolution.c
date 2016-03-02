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
double result;

double farmer(int);
void worker(int);
double generateTask(stack *tasks, double a, double b, double fa, double fb, double abarea, int numprocs);

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
    for (i = 0; i < numprocs; i++) {
      tasks_per_process[i] = 0;
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
    int i, tag, who;
    MPI_Status status;
    double temp[5];
    double *task;
    
    result = 0;
    
    generateTask(tasks, A, B, F(A), F(B), (F(A)+F(B)) * (B-A)/2, numprocs);
    
    return result;
}

double generateTask(stack *tasks, double a, double b, double fa, double fb, double abarea, int numprocs){
    double points[5], temp[5];
    double left, mid, fleft, larea;
    double right, fmid, fright, rarea;
    int i, tag, who;
    MPI_Status status;
    double *task;
    
    points[0] = a;
    points[1] = b;
    points[2] = fa;
    points[3] = fb;
    points[4] = abarea;
    push(points, tasks);
    
    printf("Push %f, %f, %f, %f, %f \n", a, b, fa, fb, abarea);
    
    i = (rand() % (numprocs-1)) + 1;
    task = pop(tasks);
    MPI_Send(task, 5, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    
    MPI_Recv(&temp, 5, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    who = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    larea = temp[0];
    rarea = temp[1];
    tasks_per_process[tag] += 1;
    
    // Create more tasks or save result
    if (temp[2] != -1 && temp[3] != -1 && temp[4] != -1){
        
        fleft = F(temp[2]);
        fmid = F(temp[3]);
        fright = F(temp[4]);
        
        generateTask(tasks, left, mid, fleft, fmid, larea, numprocs);
        generateTask(tasks, mid, right, fmid, fright, rarea, numprocs);
    }else{
        result += larea + rarea;
    }
}

void worker(int mypid) {
    // Init counters for worker
    int tasksdone = 0;
    double workdone = 0;
    
    // Init variables to receive task
    double task[5];
    int tag;
    double result[5];
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
        larea = (fleft + fmid) * (mid - left) / 2;
        rarea = (fmid + fright) * (right - mid) / 2;
        
        result[0] = larea;
        result[1] = rarea;
        
        if( fabs((larea + rarea) - lrarea) > EPSILON ) {
            result[2] = left;
            result[3] = mid;
            result[4] = right;
        }else{
            result[2] = -1;
            result[3] = -1;
            result[4] = -1;
        }
        
        // Update counters
        workdone+= result[0];
        tasksdone++;
        
        // Send result
        MPI_Send(result, 5, MPI_DOUBLE, 0, mypid, MPI_COMM_WORLD);
    
        // Receive next task
        MPI_Recv(&task, 5, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        tag = status.MPI_TAG;
    }
    printf("Worker %d solved %d tasks totalling %f units of work \n", mypid, tasksdone, workdone);
}

