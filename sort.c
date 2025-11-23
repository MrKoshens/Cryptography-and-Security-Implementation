#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define NUM_SORTS 4
#define NUM_RUNS 10000
#define MIN_SIZE 100
#define MAX_SIZE 1000
#define STEP_SIZE 100

// Holds raw metrics for a single run
typedef struct {
    unsigned long swaps;
    unsigned long comps;
    clock_t cycles;
} SortMetrics;

// Summarizes min, max, sum, avg, median
typedef struct {
    unsigned long min;
    unsigned long max;
    unsigned long long sum;
    double avg;
    double median;
    unsigned long *all_runs;  // array of length NUM_RUNS
} StatSummary;

// Comparator for qsort on unsigned long
static int compare_ulong(const void *a, const void *b) {
    unsigned long va = *(const unsigned long*)a;
    unsigned long vb = *(const unsigned long*)b;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

// Compute min, max, sum, avg, median for 'data' of length NUM_RUNS
void compute_statistics(StatSummary *stats, unsigned long *data) {
    stats->min = data[0];
    stats->max = data[0];
    stats->sum = 0;
    for (int i = 0; i < NUM_RUNS; ++i) {
        if (data[i] < stats->min) stats->min = data[i];
        if (data[i] > stats->max) stats->max = data[i];
        stats->sum += data[i];
    }
    stats->avg = stats->sum / (double)NUM_RUNS;

    // Copy into all_runs and sort for median
    memcpy(stats->all_runs, data, NUM_RUNS * sizeof(unsigned long));
    qsort(stats->all_runs, NUM_RUNS, sizeof(unsigned long), compare_ulong);
    if (NUM_RUNS % 2 == 0) {
        stats->median = (stats->all_runs[NUM_RUNS/2 - 1] + stats->all_runs[NUM_RUNS/2]) / 2.0;
    } else {
        stats->median = stats->all_runs[NUM_RUNS/2];
    }
}

// Swap helper
static void do_swap(int *a, int *b, unsigned long *swap_count) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
    (*swap_count)++;
}

void bubble_sort(int *arr, int n, SortMetrics *m) {
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            m->comps++;
            if (arr[j] > arr[j + 1])
                do_swap(&arr[j], &arr[j + 1], &m->swaps);
        }
    }
}

void merge(int *arr, int l, int m, int r, SortMetrics *met) {
    int n1 = m - l + 1;
    int n2 = r - m;
    int *L = malloc(n1 * sizeof(int));
    int *R = malloc(n2 * sizeof(int));
    memcpy(L, arr + l, n1 * sizeof(int));
    memcpy(R, arr + m + 1, n2 * sizeof(int));
    int i = 0, j = 0, k = l;
    while (i < n1 && j < n2) {
        met->comps++;
        if (L[i] <= R[j]) arr[k++] = L[i++];
        else arr[k++] = R[j++];
        met->swaps++;
    }
    while (i < n1) { arr[k++] = L[i++]; met->swaps++; }
    while (j < n2) { arr[k++] = R[j++]; met->swaps++; }
    free(L);
    free(R);
}

void merge_sort(int *arr, int l, int r, SortMetrics *m) {
    if (l < r) {
        int mid = l + (r - l) / 2;
        merge_sort(arr, l, mid, m);
        merge_sort(arr, mid + 1, r, m);
        merge(arr, l, mid, r, m);
    }
}

int partition(int *arr, int low, int high, SortMetrics *m) {
    int pivot = arr[high];
    int i = low - 1;
    for (int j = low; j < high; ++j) {
        m->comps++;
        if (arr[j] < pivot) {
            i++;
            do_swap(&arr[i], &arr[j], &m->swaps);
        }
    }
    do_swap(&arr[i+1], &arr[high], &m->swaps);
    return i + 1;
}

void quick_sort(int *arr, int low, int high, SortMetrics *m) {
    if (low < high) {
        int pi = partition(arr, low, high, m);
        quick_sort(arr, low, pi - 1, m);
        quick_sort(arr, pi + 1, high, m);
    }
}

void heapify(int *arr, int n, int i, SortMetrics *m) {
    int largest = i;
    int l = 2 * i + 1;
    int r = 2 * i + 2;
    if (l < n) { m->comps++; if (arr[l] > arr[largest]) largest = l; }
    if (r < n) { m->comps++; if (arr[r] > arr[largest]) largest = r; }
    if (largest != i) {
        do_swap(&arr[i], &arr[largest], &m->swaps);
        heapify(arr, n, largest, m);
    }
}

void heap_sort(int *arr, int n, SortMetrics *m) {
    for (int i = n/2 - 1; i >= 0; --i)
        heapify(arr, n, i, m);
    for (int i = n - 1; i > 0; --i) {
        do_swap(&arr[0], &arr[i], &m->swaps);
        heapify(arr, i, 0, m);
    }
}

void run_experiment(FILE *fout) {
    srandom(time(NULL));
    int *original = malloc(MAX_SIZE * sizeof(int));
    int *copy     = malloc(MAX_SIZE * sizeof(int));
    const char *names[NUM_SORTS] = {"Bubble", "Merge", "Quick", "Heap"};

    for (int size = MIN_SIZE; size <= MAX_SIZE; size += STEP_SIZE) {
        fprintf(fout, "Array Size: %d\n", size);

        for (int s = 0; s < NUM_SORTS; ++s) {
            unsigned long  swaps[NUM_RUNS], comps[NUM_RUNS];
            clock_t        cycles[NUM_RUNS];
            StatSummary ss_swaps = {.all_runs = malloc(NUM_RUNS * sizeof(unsigned long))};
            StatSummary ss_comps = {.all_runs = malloc(NUM_RUNS * sizeof(unsigned long))};
            StatSummary ss_cycles = {.all_runs = malloc(NUM_RUNS * sizeof(unsigned long))};

            for (int run=0; run<NUM_RUNS; ++run) {
                for (int i=0; i<size; ++i)
                    original[i] = (int)(random() % 10000);
                memcpy(copy, original, size * sizeof(int));
                SortMetrics m = {0,0,0};
                clock_t start = clock();
                switch(s) {
                    case 0: bubble_sort(copy, size, &m); break;
                    case 1: merge_sort(copy, 0, size-1, &m); break;
                    case 2: quick_sort(copy, 0, size-1, &m); break;
                    case 3: heap_sort(copy, size, &m); break;
                }
                cycles[run] = clock() - start;
                comps[run]  = m.comps;
                swaps[run]  = m.swaps;
            }

            compute_statistics(&ss_swaps, swaps);
            compute_statistics(&ss_comps, comps);
            compute_statistics(&ss_cycles, cycles);

            fprintf(fout, "%s Sort:\n", names[s]);
            fprintf(fout, "Swaps:   Min=%lu Max=%lu Median=%.2f Avg=%.2f\n",
                    ss_swaps.min, ss_swaps.max, ss_swaps.median, ss_swaps.avg);
            fprintf(fout, "Comps:   Min=%lu Max=%lu Median=%.2f Avg=%.2f\n",
                    ss_comps.min, ss_comps.max, ss_comps.median, ss_comps.avg);
            fprintf(fout, "Cycles:  Min=%lu Max=%lu Median=%.2f Avg=%.2f\n\n",
                    (unsigned long)ss_cycles.min, (unsigned long)ss_cycles.max,
                    ss_cycles.median, ss_cycles.avg);

            free(ss_swaps.all_runs);
            free(ss_comps.all_runs);
            free(ss_cycles.all_runs);
        }
        fprintf(fout, "-------------------------------\n");
    }
    free(original);
    free(copy);
}

int main() {
    FILE *fout = fopen("results.txt", "w");
    if (!fout) { perror("results.txt"); return 1; }
    run_experiment(fout);
    fclose(fout);
    return 0;
}
