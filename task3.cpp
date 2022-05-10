#include <iostream>
#include <math.h>
#include <random>
#include <fstream>
using namespace std;

double f(double x, double y, double z);
double g(double x, double y, double z);

void printArray1D(double arr[], int arr_len){
    for (int i=0; i<arr_len; i++){
        cout << arr[i] << " ";
    }
    cout << endl << endl;
}

void printArray2D(
    double** arr, 
    int row_cnt, 
    int col_cnt){
    for (int i=0; i<row_cnt; i++){
        for (int j=0; j<col_cnt; j++){
            cout << arr[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl;
}

double l1nrm(
    double* act,
    double* apr,
    int L){

    double l1 = 0;

    for (int i=0; i<L; i++){
        l1 += abs(apr[i] - act[i]);
    }

    return l1;
}

double l2nrm(
    double* act,
    double* apr,
    int L){

    double l2 = 0;

    for (int i=0; i<L; i++){
        l2 += pow(apr[i] - act[i],2);
    }

    return sqrt(l2);
}

double l0nrm(
    double* act,
    double* apr,
    int L){

    double l_inf = 0;
    double tmp = 0;

    for (int i=0; i<L; i++){
        tmp = abs(apr[i] - act[i]);
        if (l_inf<tmp){
            l_inf = abs(apr[i] - act[i]);
        }
    }

    return l_inf;
}

double* generate_zero_vector(
    double* vec,
    int L){

    for (int i=0; i<L; i++){
        vec[i] = 0;
    }

    return vec;
}

double* extractEverySecondElement(
    double arr[], 
    int arr_len) {
   
    const int new_len = (arr_len + 1) / 2;
    double* new_arr = new double[new_len];

    int j = 0;
    for (int i = 0; i < arr_len; i += 2) {
        new_arr[j] = arr[i];
        j++;
    }

    return new_arr;

}

double* applyRungeRule(
    double coarse_mesh[], // less points
    double fine_mesh[],   // more points
    int mesh_size
    ) {

    double* runge_array = new double[mesh_size];
    int acc_ord = 1;

    for (int i = 0; i < mesh_size; i++) {
        if (i == 0 || i == mesh_size - 1) {
            acc_ord = 1;
        }
        else {
            acc_ord = 2;
        }
        runge_array[i] = (fine_mesh[i] - coarse_mesh[i]) / (pow(0.5, acc_ord) - 1);
    }

    return runge_array;

}

double* generateNodes(
    double step,
    double a,
    double b){

    int n_nodes = int((b - a) / step + 1);
    double* x = new double[n_nodes];
    x[0] = a;
    for (int i=0; i<n_nodes-1; i++){
        x[i+1] = x[i] + step;
        if (x[i] == b){
            break;
        }
    }

    return x;
}

void printToFile(
    double nodes_x[],
    double nodes_y[],
    double nodes_y1[],
    int n_nodes) {

    fstream output;
    output.open("output", ios::out);

    if (!output) {
        cout << "File was not created!" << endl;
    }
    else {
        for (int i = 0; i < n_nodes; i++) {
            if (i == n_nodes - 1) {
                output << nodes_x[i];
            }
            else {
                output << nodes_x[i] << ",";
            }
        }
        output << endl;
        for (int i = 0; i < n_nodes; i++) {
            if (i == n_nodes - 1) {
                output << nodes_y[i];
            }
            else {
                output << nodes_y[i] << ",";
            }
        }
        output << endl;
        for (int i = 0; i < n_nodes; i++) {
            if (i == n_nodes - 1) {
                output << nodes_y1[i];
            }
            else {
                output << nodes_y1[i] << ",";
            }
        }
    }

}

void appendToFile(
    double nodes_x[],
    double nodes_y[],
    double nodes_y1[],
    int n_nodes) {

    fstream output;
    output.open("output", std::ios_base::app);

    if (!output) {
        cout << "File was not created!" << endl;
    }
    else {
        output << endl;
        for (int i = 0; i < n_nodes; i++) {
            if (i == n_nodes - 1) {
                output << nodes_x[i];
            }
            else {
                output << nodes_x[i] << ",";
            }
        }
        output << endl;
        for (int i = 0; i < n_nodes; i++) {
            if (i == n_nodes - 1) {
                output << nodes_y[i];
            }
            else {
                output << nodes_y[i] << ",";
            }
        }
        output << endl;
        for (int i = 0; i < n_nodes; i++) {
            if (i == n_nodes - 1) {
                output << nodes_y1[i];
            }
            else {
                output << nodes_y1[i] << ",";
            }
        }
    }

}

void euler(
    double** X, 
    double** Y, 
    double** Y1, 
    int* n_points, 
    double* step){
    
    for (int i = 0; i < 3; i++) {
        for (int j = 1; j < n_points[i]; j++)
        {
            Y[i][j]  = Y[i][j-1] + step[i] * f(X[i][j-1],Y[i][j-1],Y1[i][j-1]);
            Y1[i][j] = Y1[i][j-1] + step[i] * g(X[i][j-1],Y[i][j-1],Y1[i][j-1]);
        }
    }

}

void eulerWithRecalc(
    double** X, 
    double** Y, 
    double** Y1, 
    int* n_points, 
    double* step){
    
    double k0;
    double l0;

    for (int i = 0; i < 3; i++) {
        for (int j = 1; j < n_points[i]; j++)
        {
            k0 = Y[i][j-1] +step[i] * f(X[i][j-1],Y[i][j-1],Y1[i][j-1]);
            l0 = Y1[i][j-1] + step[i] * g(X[i][j-1],Y[i][j-1],Y1[i][j-1]);

            Y[i][j]  = Y[i][j-1] + step[i]/2 * (f(X[i][j-1],Y[i][j-1],Y1[i][j-1])+f(X[i][j],k0,l0));
            Y1[i][j] = Y1[i][j-1] + step[i]/2 * (g(X[i][j],k0,l0)+g(X[i][j-1],Y[i][j-1],Y1[i][j-1]));
        }
    }

}

void rungeKutta2(
    double** X, 
    double** Y, 
    double** Y1, 
    int* n_points, 
    double* step){
    
    double k0, k1;
    double l0, l1;

    for (int i = 0; i < 3; i++) {
        for (int j = 1; j < n_points[i]; j++)
        {
            k0 = step[i] * f(X[i][j-1],Y[i][j-1],Y1[i][j-1]);
            l0 = step[i] * g(X[i][j-1],Y[i][j-1],Y1[i][j-1]);

            k1 = step[i] * f(X[i][j-1]+step[i]/2,Y[i][j-1]+k0/2,Y1[i][j-1]+l0/2);
            l1 = step[i] * g(X[i][j-1]+step[i]/2,Y[i][j-1]+k0/2,Y1[i][j-1]+l0/2);

            Y[i][j]  = Y[i][j-1] + k1;
            Y1[i][j] = Y1[i][j-1] + l1;
        }
    }

}

void rungeKutta4(
    double** X, 
    double** Y, 
    double** Y1, 
    int* n_points, 
    double* step){
    
    double k0, k1, k2, k3;
    double l0, l1, l2, l3;

    for (int i = 0; i < 3; i++) {
        for (int j = 1; j < n_points[i]; j++)
        {
            k0 = step[i] * f(X[i][j-1],Y[i][j-1],Y1[i][j-1]);
            l0 = step[i] * g(X[i][j-1],Y[i][j-1],Y1[i][j-1]);

            k1 = step[i] * f(X[i][j-1]+step[i]/2.0,Y[i][j-1]+k0/2.0,Y1[i][j-1]+l0/2.0);
            l1 = step[i] * g(X[i][j-1]+step[i]/2.0,Y[i][j-1]+k0/2.0,Y1[i][j-1]+l0/2.0);

            k2 = step[i] * f(X[i][j-1]+step[i]/2.0,Y[i][j-1]+k1/2.0,Y1[i][j-1]+l1/2.0);
            l2 = step[i] * g(X[i][j-1]+step[i]/2.0,Y[i][j-1]+k1/2.0,Y1[i][j-1]+l1/2.0);

            k3 = step[i] * f(X[i][j-1]+step[i],Y[i][j-1]+k2,Y1[i][j-1]+l2);
            l3 = step[i] * g(X[i][j-1]+step[i],Y[i][j-1]+k2,Y1[i][j-1]+l2);

            Y[i][j]  = Y[i][j-1] + 1.0/6.0*(k0+2.0*k1+2.0*k2+k3);
            Y1[i][j] = Y1[i][j-1] + 1.0/6.0*(l0+2.0*l1+2.0*l2+l3);
        }
    }

}

void adams3(
    double** X, 
    double** Y, 
    double** Y1, 
    int* n_points, 
    double* step){
    
    double k0, k1, k2, k3;
    double l0, l1, l2, l3;

    for (int i = 0; i < 3; i++) {
        for (int j = 1; j < n_points[i]; j++)
        {
            if (j < 3){ // runge-kutta4 here
                k0 = step[i] * f(X[i][j-1],Y[i][j-1],Y1[i][j-1]);
                l0 = step[i] * g(X[i][j-1],Y[i][j-1],Y1[i][j-1]);

                k1 = step[i] * f(X[i][j-1]+step[i]/2.0,Y[i][j-1]+k0/2.0,Y1[i][j-1]+l0/2.0);
                l1 = step[i] * g(X[i][j-1]+step[i]/2.0,Y[i][j-1]+k0/2.0,Y1[i][j-1]+l0/2.0);

                k2 = step[i] * f(X[i][j-1]+step[i]/2.0,Y[i][j-1]+k1/2.0,Y1[i][j-1]+l1/2.0);
                l2 = step[i] * g(X[i][j-1]+step[i]/2.0,Y[i][j-1]+k1/2.0,Y1[i][j-1]+l1/2.0);

                k3 = step[i] * f(X[i][j-1]+step[i],Y[i][j-1]+k2,Y1[i][j-1]+l2);
                l3 = step[i] * g(X[i][j-1]+step[i],Y[i][j-1]+k2,Y1[i][j-1]+l2);

                Y[i][j]  = Y[i][j-1] + 1.0/6.0*(k0+2.0*k1+2.0*k2+k3);
                Y1[i][j] = Y1[i][j-1] + 1.0/6.0*(l0+2.0*l1+2.0*l2+l3);
            }
            else{
                k0 = step[i] * f(X[i][j-1],Y[i][j-1],Y1[i][j-1]);
                l0 = step[i] * g(X[i][j-1],Y[i][j-1],Y1[i][j-1]);

                k1 = step[i] * f(X[i][j-2],Y[i][j-2],Y1[i][j-2]);
                l1 = step[i] * g(X[i][j-2],Y[i][j-2],Y1[i][j-2]);

                k2 = step[i] * f(X[i][j-3],Y[i][j-3],Y1[i][j-3]);
                l2 = step[i] * g(X[i][j-3],Y[i][j-3],Y1[i][j-3]);

                Y[i][j]  = Y[i][j-1] + 1.0/12.0*(23.0*k0-16.0*k1+5.0*k2);
                Y1[i][j] = Y1[i][j-1] + 1.0/12.0*(23.0*l0-16.0*l1+5.0*l2);
            }
        }
    }

}

void rungeRelative(
    double** rungeMatrix,
    double** Y,
    double** zeros,
    int* n_points,
    int method_n
){

    double* runge_array = new double [n_points[0]];
    double runge_l1, runge_l2, runge_l0;
    double runge_l1_rel, runge_l2_rel, runge_l0_rel;

    runge_array  = applyRungeRule(Y[0], extractEverySecondElement(Y[1], n_points[1]), n_points[0]);
    runge_l1     = l1nrm(runge_array, zeros[0], n_points[0]);
    runge_l2     = l2nrm(runge_array, zeros[0], n_points[0]);
    runge_l0     = l0nrm(runge_array, zeros[0], n_points[0]);
    runge_l1_rel = runge_l1 / l1nrm(Y[0], zeros[0], n_points[0]);
    runge_l2_rel = runge_l2 / l2nrm(Y[0], zeros[0], n_points[0]);
    runge_l0_rel = runge_l0 / l0nrm(Y[0], zeros[0], n_points[0]);

    rungeMatrix[method_n-1][0]=runge_l0_rel;
    rungeMatrix[method_n-1][1]=runge_l1_rel;
    rungeMatrix[method_n-1][2]=runge_l2_rel;
}

void deviations(
    double** deviationsMatrix,
    double** Y,
    double** zeros,
    int* n_points,
    int method_n
){

    double abs_l1_1 = 0, abs_l1_2 = 0;
    double abs_l2_1 = 0, abs_l2_2 = 0;
    double abs_l0_1 = 0, abs_l0_2 = 0;

    double rel_l1_1 = 0, rel_l1_2 = 0;
    double rel_l2_1 = 0, rel_l2_2 = 0;
    double rel_l0_1 = 0, rel_l0_2 = 0;

    abs_l0_1 = l0nrm(Y[0], extractEverySecondElement(Y[1], n_points[1]), n_points[0]);
    abs_l1_1 = l1nrm(Y[0], extractEverySecondElement(Y[1], n_points[1]), n_points[0]);
    abs_l2_1 = l2nrm(Y[0], extractEverySecondElement(Y[1], n_points[1]), n_points[0]);

    rel_l0_1 = abs_l0_1 / l0nrm(Y[0], zeros[0], n_points[0]);
    rel_l1_1 = abs_l1_1 / l1nrm(Y[0], zeros[0], n_points[0]);
    rel_l2_1 = abs_l2_1 / l2nrm(Y[0], zeros[0], n_points[0]);

    abs_l0_2 = l0nrm(extractEverySecondElement(Y[0], n_points[0]), Y[2], n_points[2]);
    abs_l1_2 = l1nrm(extractEverySecondElement(Y[0], n_points[0]), Y[2], n_points[2]);
    abs_l2_2 = l2nrm(extractEverySecondElement(Y[0], n_points[0]), Y[2], n_points[2]);

    rel_l0_2 = abs_l0_2 / l0nrm(extractEverySecondElement(Y[0], n_points[0]), zeros[2], n_points[2]);
    rel_l1_2 = abs_l1_2 / l1nrm(extractEverySecondElement(Y[0], n_points[0]), zeros[2], n_points[2]);
    rel_l2_2 = abs_l2_2 / l2nrm(extractEverySecondElement(Y[0], n_points[0]), zeros[2], n_points[2]);

    deviationsMatrix[6*(method_n-1)+0][0]=abs_l0_1;
    deviationsMatrix[6*(method_n-1)+1][0]=abs_l1_1;
    deviationsMatrix[6*(method_n-1)+2][0]=abs_l2_1;
    deviationsMatrix[6*(method_n-1)+3][0]=rel_l0_1;
    deviationsMatrix[6*(method_n-1)+4][0]=rel_l1_1;
    deviationsMatrix[6*(method_n-1)+5][0]=rel_l2_1;

    deviationsMatrix[6*(method_n-1)+0][1]=abs_l0_2;
    deviationsMatrix[6*(method_n-1)+1][1]=abs_l1_2;
    deviationsMatrix[6*(method_n-1)+2][1]=abs_l2_2;
    deviationsMatrix[6*(method_n-1)+3][1]=rel_l0_2;
    deviationsMatrix[6*(method_n-1)+4][1]=rel_l1_2;
    deviationsMatrix[6*(method_n-1)+5][1]=rel_l2_2;
}

void printRunge(
    double** rungeMatrix,
    double** rungeMatrix1
){

    cout << "\n \t\t---> Relative Runge Error for y(x) <---\n\n";
    cout << "       \t\t\t|\t l0 \t|\t l1 \t|\t l2\n";
    cout << "----------------------------------------------------------------------\n";

    for (int i = 1; i < 6; i++){
        
        switch(i){
            case 1 : cout << "Euler \t\t\t|"; break;
            case 2 : cout << "Euler w/ recalc \t|"; break;
            case 3 : cout << "Runge-Kutta (2o) \t|"; break;
            case 4 : cout << "Runge-Kutta (4o) \t|"; break;
            case 5 : cout << "Adams \t\t\t|"; break;
        }

        printf(" %e  | ", rungeMatrix[i-1][0]);
        printf(" %e | ", rungeMatrix[i-1][1]);
        printf(" %e  \n", rungeMatrix[i-1][2]);
    }

    cout << "\n\n \t\t---> Relative Runge Error for y'(x) <---\n\n";
    cout << "       \t\t\t|\t l0 \t|\t l1 \t|\t l2\n";
    cout << "----------------------------------------------------------------------\n";

    for (int i = 1; i < 6; i++){
        
        switch(i){
            case 1 : cout << "Euler \t\t\t|"; break;
            case 2 : cout << "Euler w/ recalc \t|"; break;
            case 3 : cout << "Runge-Kutta (2o) \t|"; break;
            case 4 : cout << "Runge-Kutta (4o) \t|"; break;
            case 5 : cout << "Adams \t\t\t|"; break;
        }

        printf(" %e  | ", rungeMatrix1[i-1][0]);
        printf(" %e | ", rungeMatrix1[i-1][1]);
        printf(" %e  \n", rungeMatrix1[i-1][2]);
    }
    cout << endl;
}

void printDeviations(
    double** deviationsMatrix,
    double** deviationsMatrix1
){

    cout << "\n \t\t   ---> DEVIATION of y(x) <---\n"     << endl; 
    cout << "                 \t\t |\t  H/2 \t    |\t    2*H\n";
    for (int i = 1; i < 6; i++){     
        cout << "------------------------------------------------------------------------\n";   
        switch(i){
            case 1 : cout << "Euler \t\t\t|"; break;
            case 2 : cout << "Euler w/ recalc \t|"; break;
            case 3 : cout << "Runge-Kutta (2o) \t|"; break;
            case 4 : cout << "Runge-Kutta (4o) \t|"; break;
            case 5 : cout << "Adams \t\t\t|"; break;
        }
        for (int j = 0; j < 3; j++){
            if (j==0){
                printf(" l%i_abs", j);
            }
            else{
                printf("      \t\t\t| l%i_abs", j);
            }
            printf(" |   %e   |   %e\n", deviationsMatrix[6*(i-1)+j][0], deviationsMatrix[6*(i-1)+j][1]);
        }
        cout << "\t\t\t------------------------------------------------\n";
        for (int j = 3; j < 6; j++){
            printf("      \t\t\t| l%i_rel", j-3);
            printf(" |   %e   |   %e\n", deviationsMatrix[6*(i-1)+j][0], deviationsMatrix[6*(i-1)+j][1]);
        }
    }
    
    cout << endl;

    cout << "\n \t\t   ---> DEVIATION of y'(x) <---\n"     << endl; 
    cout << "                 \t\t |\t  H/2 \t    |\t    2*H\n";
    for (int i = 1; i < 6; i++){     
        cout << "------------------------------------------------------------------------\n";   
        switch(i){
            case 1 : cout << "Euler \t\t\t|"; break;
            case 2 : cout << "Euler w/ recalc \t|"; break;
            case 3 : cout << "Runge-Kutta (2o) \t|"; break;
            case 4 : cout << "Runge-Kutta (4o) \t|"; break;
            case 5 : cout << "Adams \t\t\t|"; break;
        }
        for (int j = 0; j < 3; j++){
            if (j==0){
                printf(" l%i_abs", j);
            }
            else{
                printf("      \t\t\t| l%i_abs", j);
            }
            printf(" |   %e   |   %e\n", deviationsMatrix1[6*(i-1)+j][0], deviationsMatrix1[6*(i-1)+j][1]);
        }
        cout << "\t\t\t------------------------------------------------\n";
        for (int j = 3; j < 6; j++){
            printf("      \t\t\t| l%i_rel", j-3);
            printf(" |   %e   |   %e\n", deviationsMatrix1[6*(i-1)+j][0], deviationsMatrix1[6*(i-1)+j][1]);
        }
    }

    cout << endl;
}


double f(double x, double y, double z){ // dy/dx = f(x,y,z)
    return z;
}

double g(double x, double y, double z){ // dz/dx = g(x,y,z)
    return z*exp(x+2)+y*sin(x)+tan(x+1);
    //return -sin(x)*x*cos(x);
}


int main(){

    double a  = 0;
    double b  = 1;
    double y0 = -5; // y(0)
    double z0 = 5;  // y'(0)

    int n_intervals = 100;

    double* step        = new double[3];
    int* n_points       = new int[3];
    double** X          = new double*[3];
    double** Y          = new double*[3];
    double** Y1         = new double*[3];
    double** zeros      = new double*[3];

    double** rungeMatrix  = new double*[5];  // method -> metric(relative)
    double** rungeMatrix1 = new double*[5];

    double** deviationsMatrix  = new double*[30];  // method -> metric(mesh)
    double** deviationsMatrix1 = new double*[30];

    step[0] = (b - a) / n_intervals; // initial mesh
    step[1] = step[0]/2;             // mesh with more points
    step[2] = step[0]*2;             // mesh with less points

    // calculate number of points
    for (int i = 0; i < 3; i++) {
        n_points[i] = int((b - a) / step[i] + 1);
    }

    //generate meshes and set initial values
    for (int i = 0; i < 3; i++) {
        X[i]     = new double[n_points[i]];
        Y[i]     = new double[n_points[i]];
        Y1[i]    = new double[n_points[i]];
        zeros[i] = new double[n_points[i]];

        X[i]     = generateNodes(step[i], a, b);
        zeros[i] = generate_zero_vector(zeros[i], n_points[i]);
        Y[i][0]  = y0;
        Y1[i][0] = z0;
    }

    for (int i = 0; i < 5; i++) {
        rungeMatrix[i]  = new double[3];
        rungeMatrix1[i] = new double[3];
    }

    for (int i = 0; i < 30; i++) {
        deviationsMatrix[i]  = new double[2];
        deviationsMatrix1[i] = new double[2];
    }

    // Euler method
    euler(X, Y, Y1, n_points, step);
    printToFile(X[0], Y[0], Y1[0], n_points[0]);
    rungeRelative(rungeMatrix, Y, zeros, n_points, 1);
    rungeRelative(rungeMatrix1, Y1, zeros, n_points, 1);
    deviations(deviationsMatrix, Y, zeros, n_points, 1);
    deviations(deviationsMatrix1, Y1, zeros, n_points, 1);

    // Euler method with recalculation
    eulerWithRecalc(X, Y, Y1, n_points, step);
    appendToFile(X[0], Y[0], Y1[0], n_points[0]);
    rungeRelative(rungeMatrix, Y, zeros, n_points, 2);
    rungeRelative(rungeMatrix1, Y1, zeros, n_points, 2);
    deviations(deviationsMatrix, Y, zeros, n_points, 2);
    deviations(deviationsMatrix1, Y1, zeros, n_points, 2);

    // Runge-Kutta method (2nd order)
    rungeKutta2(X, Y, Y1, n_points, step);
    appendToFile(X[0], Y[0], Y1[0], n_points[0]);
    rungeRelative(rungeMatrix, Y, zeros, n_points, 3);
    rungeRelative(rungeMatrix1, Y1, zeros, n_points, 3);
    deviations(deviationsMatrix, Y, zeros, n_points, 3);
    deviations(deviationsMatrix1, Y1, zeros, n_points, 3);

    // Runge-Kutta method (4th order)
    rungeKutta4(X, Y, Y1, n_points, step);
    appendToFile(X[0], Y[0], Y1[0], n_points[0]);
    rungeRelative(rungeMatrix, Y, zeros, n_points, 4);
    rungeRelative(rungeMatrix1, Y1, zeros, n_points, 4);
    deviations(deviationsMatrix, Y, zeros, n_points, 4);
    deviations(deviationsMatrix1, Y1, zeros, n_points, 4);

    // Adams method (3rd order)
    adams3(X, Y, Y1, n_points, step);
    appendToFile(X[0], Y[0], Y1[0], n_points[0]);
    rungeRelative(rungeMatrix, Y, zeros, n_points, 5);
    rungeRelative(rungeMatrix1, Y1, zeros, n_points, 5);
    deviations(deviationsMatrix, Y, zeros, n_points, 5);
    deviations(deviationsMatrix1, Y1, zeros, n_points, 5);

    // Deviations from mesh H
    printDeviations(deviationsMatrix, deviationsMatrix1);

    // Runge's estimation
    printRunge(rungeMatrix, rungeMatrix1);

    return 0;
}