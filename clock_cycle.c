#include <stdio.h>
#include <stdint.h>
#include <mach/mach_time.h>

// Euclidean GCD
unsigned int gcd(unsigned int a, unsigned int b) {
    while (b != 0) {
        unsigned int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

int main() {
    unsigned int a = 2527364, b = 1058273;

    // Get timebase info for conversion
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);

    uint64_t start = mach_absolute_time();

    unsigned int result = gcd(a, b);

    uint64_t end = mach_absolute_time();

    // Convert to nanoseconds
    uint64_t elapsed_ns = (end - start) * timebase.numer / timebase.denom;

    // Estimated cycles: M1 Performance core ~3.2 GHz
    double estimated_cycles = (double)elapsed_ns * 3.2;

    printf("GCD(%u, %u) = %u\n", a, b, result);
    printf("Time taken: %llu ns\n", elapsed_ns);
    printf("Estimated clock cycles: %.0f\n", estimated_cycles);

    return 0;
}