/**
 * PPLS Assignment 2
 * Maria Velarde (s1556573)
 *
 * TO COMPILE
 * /usr/lib64/openmpi/bin/mpicc -o aquadSolution aquadSolution.c -lm
 *
 * TO RUN
 *  /usr/lib64/openmpi/bin/mpirun -c 5 ./aquadSolution
 * 
 * Implementation Strategy
 * -----------------------
 * Using MPI, the recursive algorithm is simulated.
 * Given the 5 parameters A, B, F(A), F(B) and ABarea the farmer
 * will add tasks to the bag. 
 * The farmer receives from the workers either larea, rarea, left, mid, 
 * right or only larea, rarea. If five parameters are received, then a new
 * task must be added to the bag, otherwise larea and rarea are 
 * added to the final result.
 *
 * Given 5 parameters A,B, F(A), F(B) and ABarea, each
 * worker must compute mid point and left and right areas.
 * While fabs((larea + rarea) - ABarea) is bigger than a
 * threshold, the worker will return 5 parameters, otherwise
 * will only return larea and rarea.
 *
 * MPI Primitives used
 * -------------------
 * The following primitives were used because the were necessary
 * for the MPI environment to work: initialize, finalize, 
 * get number of processes and get rank of current process.
 *
 * int MPI_Init(int *argc, char ***argv)
 * -------------------------------------
 * Initializes the MPI environment.
 *
 * int MPI_Comm_size(MPI_Comm comm, int *size)
 * -------------------------------------------
 * Given a communicator, MPI_COMM_WORLD in this case,
 * gets the number of processes available.
 *
 * int MPI_Comm_rank(MPI_Comm comm, int *rank)
 * -------------------------------------------
 * Gets the rank of the current process in the communicator's group.
 *
 * int MPI_Finalize()
 * ------------------
 * Ends and clean the MPI environment
 *
 *--------------------------------------------------------------------------------
 * The following primitives were used because they suited best the implementation.
 *
 * FARMER
 * ------
 * int MPI_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
 * int dest, int sendtag, void *recvbuf, int recvcount,
 * MPI_Datatype recvtype, int source, int recvtag,
 * MPI_Comm comm, MPI_Status *status)
 * ---------------------------------------------------------------------
 * To be able to finalize the MPI environment every send must have its corresponding
 * receive. 
 * In the while, whenever a task is popped from the stack it is send to a
 * random process and the while waits for its result. 
 * Instead of using MPI_Send followed by its correspondent MPI_Recv, 
 * MPI_Sendrecv is used.
 *
 *
 * int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest,
 * int tag, MPI_Comm comm)
 * -------------------------------------------------------------------
 * Later on, the farmer needs to tell every worker that there are no
 * more tasks. MPI_Send is sufficient for this case, using MPI_Isend 
 * for example would only complicate things by creating the need to
 * probe communication. Non-blocking communication shouldn't be used
 * unless is necessary.
 *
 *
 * WORKER
 * ------
 * The worker wil receive from the farmer either many MPI_Sendrecv
 * with data or a single MPI_Send with the NO MORE TASKS tag.
 *
 * int MPI_Recv(void *buf, int count, MPI_Datatype datatype,
 * int source, int tag, MPI_Comm comm, MPI_Status *status)
 * ---------------------------------------------------------
 * The first receive corresponds to the first task received from the 
 * farmer's MPI_Sendrecv. Once in the while the receive will correspond
 * either to the farmer's MPI_Sendrecv or to the farmer's MPI_Send.
 * Depending on the tag received the worker will stay or not in the while.
 *
 * int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest,
 * int tag, MPI_Comm comm)
 * -------------------------------------------------------------------
 * After every worker's MPI_Recv corresponding to a MPI_Sendrecv from the 
 * farmer, results should be sent back.
 *
 * For readibility and not associate the MPI_Send corresponding to the results
 * of one task, with the MPI_Recv corresponding to the next task (or termination 
 * condition), separate MPI_Recv and MPI_Send were used in the worker.
 *
 *
 * References
 * ----------
 * Course overheads (http://www.inf.ed.ac.uk/teaching/courses/ppls/pplsslides.pdf)
 * MPI Documentation (https://www.open-mpi.org/doc/v1.4/)
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>


#define EPSILON 1e-3
#define F(arg)  cosh(arg)*cosh(arg)*cosh(arg)*cosh(arg)
#define A 0.0
#define B 5.0

#define SLEEPTIME 1

#define NO_MORE_TASKS 100

int *tasks_per_process;

double farmer(int);
void worker(int);

/* STACK DECLARATIONS */
typedef struct stack_node_tag stack_node;
typedef struct stack_tag stack;

struct stack_node_tag { double data[5]; stack_node *next;};
struct stack_tag { stack_node *top; };

/* stack init / destroy */
stack *new_stack();
void free_stack(stack *);

/* stack methods */
void push(double *, stack *);
double *pop (stack *);
int is_empty (stack *);

/* INT STACK DECLARATIONS */
typedef struct stack_node_tag_int stack_node_int;
typedef struct stack_tag_int stack_int;

struct stack_node_tag_int { int data; stack_node_int *next;};
struct stack_tag_int { stack_node_int *top; };

/* stack init / destroy */
stack_int *new_stack_int();
void free_stack_int(stack_int *);

/* stack methods */
void push_int(int, stack_int *);
int pop_int (stack_int *);
int is_empty_int (stack_int *);


/* Bag of tasks implementation */
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
    stack_int *idleWorkers = new_stack_int();
    
    // Init variables to receive info from MPI_Recv
    int i, tag, who, idleCount;
    MPI_Status status;
    double result;
    double points[5],temp[5];
    double *task;
    
    double left, right, fleft, fright;
    double mid, fmid, larea, rarea;
    
    // Init idle workers stacks
    idleCount = numprocs - 1;
    for (i = 1; i < numprocs; i++) {
        push_int(i, idleWorkers);
    }
    
    // Generate first task
    points[0] = A;
    points[1] = B;
    points[2] = F(A);
    points[3] = F(B);
    points[4] = (F(A)+F(B)) * (B-A)/2;
    push(points, tasks);
    
    // While there are tasks to process
    while (!is_empty(tasks)) {

        // If tasks to do and idle workers, send tasks
        while (!is_empty_int(idleWorkers) && !is_empty(tasks)) {
            int worker = pop_int(idleWorkers);
            
            // Remove from the stack
            task = pop(tasks);
            
            // Send task
            //printf("Send %f, %f, %f, %f, %f to %d \n", task[0], task[1], task[2], task[3], task[4], worker);
            // Args sent: task buffer, size of send buffer, data type, destination, origin (tag), common world
            MPI_Send(task, 5, MPI_DOUBLE, worker, 0, MPI_COMM_WORLD);
            idleCount--;
        }

        // While busy workers, receive results
        while ((numprocs - 1) > idleCount){
            // Args sent: result buffer, size of result buffer, data type, source, tag, common world, status
            MPI_Recv(&temp, 5, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            who = status.MPI_SOURCE;

            // Get data from temp
            larea = temp[0];
            rarea = temp[1];
            left = temp[2];
            mid = temp[3];
            right = temp[4];
            
            // Update amount of tasks processed by thread
            tasks_per_process[who] += 1;
            
            // Update idle workers
            push_int(who, idleWorkers);
            idleCount++;
            
            // Create more tasks or save result
            if (left != -1 && mid != -1 && right != -1){
                // Generate values for next two iterations and store in stack
                fleft = F(temp[2]);
                fmid = F(temp[3]);
                fright = F(temp[4]);
                
                points[0] = left;
                points[1] = mid;
                points[2] = fleft;
                points[3] = fmid;
                points[4] = larea;
                push(points, tasks);
                
                //printf("Farmer pushes %f, %f, %f, %f, %f into stack \n", points[0], points[1], points[2], points[3], points[4]);
                
                points[0] = mid;
                points[1] = right;
                points[2] = fmid;
                points[3] = fright;
                points[4] = rarea;
                push(points, tasks);
                //printf("Farmer pushes %f, %f, %f, %f, %f into stack \n", points[0], points[1], points[2], points[3], points[4]);
            }else{
                //printf("Farmer adds %f, %f to %f \n", larea, rarea, result);
                result += larea + rarea;
            }
        }
    }
    
    // Tell workers there are no more tasks to process
    for (i=0; i < (numprocs-1); i++) {
        MPI_Send(task, 5, MPI_DOUBLE, i+1, NO_MORE_TASKS, MPI_COMM_WORLD);
    }
    
    return result;
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
        usleep(SLEEPTIME);
        
        // Get variables
        left = task[0];
        right = task[1];
        fleft = task[2];
        fright = task[3];
        lrarea = task[4];
        
        mid = (left + right) / 2.0;
        fmid = F(mid);
        
        // Get result
        larea = (fleft + fmid) * (mid - left) / 2.0;
        rarea = (fmid + fright) * (right - mid) / 2.0;
        
        result[0] = larea;
        result[1] = rarea;
        
        //printf("Worker %d receives %f, %f, %f, %f, %f and computes %f, %f \n", mypid, left, right, fleft, fright, lrarea, larea, rarea);
        
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
        workdone+= result[0] + result[1];
        tasksdone++;
        
        // Send result
        MPI_Send(result, 5, MPI_DOUBLE, 0, mypid, MPI_COMM_WORLD);
    
        // Receive next task
        MPI_Recv(&task, 5, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        tag = status.MPI_TAG;
    }
    //printf("Worker %d solved %d tasks totalling %f units of work \n", mypid, tasksdone, workdone);
}

/* STACK IMPLEMENTATION */
// creating a new stack
stack * new_stack()
{
    stack *n;
    
    n = (stack *) malloc (sizeof(stack));
    
    n->top = NULL;
    
    return n;
}

// cleaning up after use
void free_stack(stack *s)
{
    free(s);
}

// Push data to stack s, data has to be an array of 2 doubles
void push (double *data, stack *s)
{
    stack_node *n;
    n = (stack_node *) malloc (sizeof(stack_node));
    n->data[0]  = data[0];
    n->data[1]  = data[1];
    n->data[2]  = data[2];
    n->data[3]  = data[3];
    n->data[4]  = data[4];
    
    if (s->top == NULL) {
        n->next = NULL;
        s->top  = n;
    } else {
        n->next = s->top;
        s->top = n;
    }
}

// Pop data from stack s
double * pop (stack * s)
{
    stack_node * n;
    double *data;
    
    if (s == NULL || s->top == NULL) {
        return NULL;
    }
    n = s->top;
    s->top = s->top->next;
    data = (double *) malloc(5*(sizeof(double)));
    data[0] = n->data[0];
    data[1] = n->data[1];
    data[2] = n->data[2];
    data[3] = n->data[3];
    data[4] = n->data[4];
    
    free (n);
    
    return data;
}

// Check for an empty stack
int is_empty (stack * s) {
    return (s == NULL || s->top == NULL);
}

/* INT STACK IMPLEMENTATION */
// creating a new stack
stack_int * new_stack_int()
{
    stack_int *n;
    
    n = (stack_int *) malloc (sizeof(stack_int));
    
    n->top = NULL;
    
    return n;
}

// cleaning up after use
void free_stack_int(stack_int *s)
{
    free(s);
}

// Push data to stack s, data has to be an array of 2 doubles
void push_int (int data, stack_int *s)
{
    stack_node_int *n;
    n = (stack_node_int *) malloc (sizeof(stack_node_int));
    n->data  = data;
    
    if (s->top == NULL) {
        n->next = NULL;
        s->top  = n;
    } else {
        n->next = s->top;
        s->top = n;
    }
}

// Pop data from stack s
int pop_int (stack_int * s)
{
    stack_node_int * n;
    int data;
    
    if (s == NULL || s->top == NULL) {
        return -1;
    }
    n = s->top;
    s->top = s->top->next;
    data = n->data;
    
    free (n);
    
    return data;
}

// Check for an empty stack
int is_empty_int (stack_int * s) {
    return (s == NULL || s->top == NULL);
}

