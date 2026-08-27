#define CONTROL

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>



double integrate_simpson38(int lower, int upper, double* f) {
    double integration = 0;
    integration = f[lower] + f[upper];
    for (int i = 1; i < (upper - lower); i++) {
        int k = lower + i;

        if (i % 3 == 0) {
            integration += 2 * f[k];
        }
        else {
            integration += 3 * f[k];
        }

    }

    integration *= 3.0 / 8.0;
    return integration;
}

double integrate_simpson38_with_reminder(int lower, int upper, double* f) {

    int reminder = (upper - lower) % 3;
    int final_38 = (upper - lower) - reminder;

    double integration = f[lower] + f[final_38];
    if (upper - lower % 3)
        for (int i = 1; i < final_38; i++) {
            int k = lower + i;

            if (i % 3 == 0) {
                integration += 2 * f[k];
            }
            else {
                integration += 3 * f[k];
            }

        }
    integration *= 3.0 / 8.0;

    double sum = 0.0;
    if (reminder != 0) {
        sum = f[final_38] + f[upper];
        for (int i = 1; i < upper - final_38; ++i) {
            int k = final_38 + i;
            sum += 2.0 * f[k];
        }
        sum /= 2.0;
    }


    return integration + sum;
}

// Helper function to fill our array with f(x) = x^3 values
void fill_function_values(double* f, int n, double lower_x, double h) {
    for (int i = 0; i <= n; i++) {
        double x = lower_x + i * h;
        f[i] = x * x * x * x; // f(x) = x^3
    }
}

int main() {
    // double exact_answer = 20.25; // Analytical integral of x^3 from 0 to 3
    double exact_answer = pow(3, 5) / 5.0; // Analytical integral of x^3 from 0 to 3
    double a = 0.0;              // Lower integration bound
    double b = 3.0;              // Upper integration bound

    printf("--- SIMPSON'S 3/8 RULE TEST PROGRAM ---\n");
    printf("Integrating f(x) = x^3 from %.1f to %.1f (Exact Answer = %.4f)\n\n", a, b, exact_answer);

    // =========================================================================
    // TEST 1: Number of intervals IS a multiple of 3 (n = 6)
    // =========================================================================
    int n1 = 12;
    double h1 = (b - a) / n1; // h = 0.5
    double* f1 = (double*)malloc((n1 + 1) * sizeof(double));

    fill_function_values(f1, n1, a, h1);

    // Multiply by h1 outside the function as per your setup
    double result1 = integrate_simpson38(0, n1, f1) * h1;
    double error1 = fabs(exact_answer - result1);

    printf("[Test 1] Number of intervals: %d (Valid multiple of 3)\n", n1);
    printf("         Calculated: %10.6f\n", result1);
    printf("         Abs Error : %10.6f\n\n", error1);

    free(f1);

    // =========================================================================
    // TEST 2: Number of intervals IS NOT a multiple of 3 (n = 5)
    // =========================================================================
    int n2 = 5;
    double h2 = (b - a) / n2; // h = 0.6
    double* f2 = (double*)malloc((n2 + 1) * sizeof(double));

    fill_function_values(f2, n2, a, h2);

    // Force the non-multiple of 3 into the function
    double result2 = integrate_simpson38(0, n2, f2) * h2;
    double error2 = fabs(exact_answer - result2);

    printf("[Test 2] Number of intervals: %d (NOT a multiple of 3)\n", n2);
    printf("         Calculated: %10.6f\n", result2);
    printf("         Abs Error : %10.6f\n", error2);
    printf("         *Note: The structural weight mismatch causes this error.\n");

    free(f2);

    for (int n = 3;n < 30;n++) {
        // =========================================================================
        // TEST 3: Number of intervals IS NOT a multiple of 3 (n = 5)
        // =========================================================================
        int n2 = n;
        double h2 = (b - a) / n2; // h = 0.6
        double* f2 = (double*)malloc((n2 + 1) * sizeof(double));

        fill_function_values(f2, n2, a, h2);

        // Force the non-multiple of 3 into the function
        double result2 = integrate_simpson38_with_reminder(0, n2, f2) * h2;
        double error2 = fabs(exact_answer - result2);

        // printf("[Test 2] Number of intervals: %d (NOT a multiple of 3)\n", n2);
        // printf("         Calculated: %10.6f\n", result2);
        printf("n=%d       Abs Error : %10.6f\n", n, error2);
        // printf("         *Note: The structural weight mismatch causes this error.\n");

        free(f2);
    }

    return 0;
}