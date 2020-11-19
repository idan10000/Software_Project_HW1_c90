#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


void algorithm(int N, int K, int d, int MAX_ITER, float obs[N][d], float centroids[K][d], int clusters[N],
               float sums[N][d]);

int main(int argc, char *argv[]) {
    const K = atoi(argv[1]);
    const N = atoi(argv[2]);
    const d = atoi(argv[3]);
    const MAX_ITER = atoi(argv[4]);

    assert(K > 0 && N > 0 && d > 0 && MAX_ITER > 0 && K < N);

    // declare array variables
    float (*obs)[d] = malloc(sizeof(float) * N * d);
    assert(obs != NULL);
    float (*centroids)[d] = malloc(sizeof(float) * K * d);
    assert(centroids != NULL);
    float (*sums)[d] = malloc(K * d * sizeof(float));
    assert(sums != NULL);
    int *clusters = malloc(N * sizeof(int));
    assert(clusters != NULL);

    algorithm(N, K, d, MAX_ITER, obs, centroids, clusters, sums);

    //Printing centroids
    int i, j;
    for (i = 0; i < K; ++i) {
        for (j = 0; j < d-1; ++j) {
            printf("%.2f%c", centroids[i][j],',');
        }
        printf("%.2f%c", centroids[i][d-1],'\n');
    }

    // free memory
    free(obs);
    free(centroids);
    free(clusters);
    free(sums);
    return 0;
}

void initObs(int N, int K, int d, float obs[N][d], float centroids[K][d]) {
    char c;
    float f;
    int i, j;
    for (i = 0; i < N; ++i) {
        for (j = 0; j < d; ++j) {
            assert(scanf("%f%c", &f, &c) == 2);
            if (i < K)
                centroids[i][j] = f;
            obs[i][j] = f;
        }
    }
}


double norm(int N, int d, float x[d], float cluster[N]) {
    double sum = 0;
    int i;
    for (i = 0; i < d; ++i) {
        sum += (x[i] - cluster[i]) * (x[i] - cluster[i]);
    }
    return sum;
}

int assignCluster(int N, int K, int d, float x[d], float centroids[K][d]) {
    double sum = norm(N, d, x, centroids[0]);
    double tempSum;
    int minCluster = 0;
    int i;
    for (i = 1; i < K; ++i) {
        tempSum = norm(N, d, x, centroids[i]);
        if (tempSum < sum) {
            sum = tempSum;
            minCluster = i;
        }
    }
    return minCluster;
}

void assignAllObservations(int N, int K, int d, float obs[N][d], float centroids[K][d], int clusters[N]) {
    int i;
    for (i = 0; i < N; ++i) {
        clusters[i] = assignCluster(N, K, d, obs[i], centroids);
    }
}

void resetSums(int N, int d, float sums[N][d]) {
    int i, j;
    for (i = 0; i < N; ++i)
        for (j = 0; j < d; ++j)
            sums[i][j]=0;
}

/**
 * after the clusters were updated, we update the centroids in a loop over the observations,
 * adding each observation coordinate to the corresponding centroids matrix coordinate.
 *
 * @param obs A matrix of the size N x d which represents the N d-dimensional vectors of the observations
 * @param centroids A matrix of the size K x d which represents the K d-dimensional vectors  of the centroids
 * @param clusters An array of the size N which represents the cluster that each observation is assigned to.
 * @param K The amount of clusters
 * @param N The amount of observations
 * @param d The size of the vector of each observation and centroid
 * @return If any of the centroids were updated
 */
char updateCentroids(int N, int K, int d, float obs[N][d], float centroids[K][d], int clusters[N], float sums[N][d]) {
    resetSums(N,d,sums);
    int clusterSizes[K];
    char changedAny = 0;
    float tempCentroid[d];
    int i, j;

    for (i = 0; i < K; ++i) {
        clusterSizes[i] = 0;
    }

    for (i = 0; i < N; ++i) {
        for (j = 0; j < d; ++j) {
            sums[clusters[i]][j] += obs[i][j];
        }
        clusterSizes[clusters[i]]++;
    }
    for (i = 0; i < K; ++i) {
        for (j = 0; j < d; ++j) {
            tempCentroid[j] = sums[i][j] / clusterSizes[i];
            if (!changedAny && tempCentroid[j] != centroids[i][j])
                changedAny = 1;
            centroids[i][j] = tempCentroid[j];
        }
    }
    return changedAny;
}


void algorithm(int N, int K, int d, int MAX_ITER, float obs[N][d], float centroids[K][d], int clusters[N],
               float sums[N][d]) {
    initObs(N, K, d, obs, centroids);
    char changedCluster = 1;
    int i;
    for (i = 0; i < MAX_ITER && changedCluster; ++i) {
        assignAllObservations(N, K, d, obs, centroids, clusters);
        changedCluster = updateCentroids(N, K, d, obs, centroids, clusters, sums);
    }
}
