#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


void algorithm(int MAX_ITER, float **obs, float **centroids, int *clusters, float **sums);


int K, N, d;

int main(int argc, char *argv[]) {
    int MAX_ITER;
    int i;
    float **obs, **centroids, **sums, **cobs, **ccentroids, **csums;
    int *clusters, *cclusters;
    assert(argc == 5);
    K = atoi(argv[1]);
    N = atoi(argv[2]);
    d = atoi(argv[3]);
    MAX_ITER = atoi(argv[4]);


    assert(K > 0 && N > 0 && d > 0 && MAX_ITER > 0 && K < N);

    /* declare array variables */
    obs = malloc(sizeof(float *) * N);
    assert(obs != NULL);
    for (i = 0; i < N; ++i) {
        obs[i] = malloc(sizeof(float) * d);
        assert(obs[i] != NULL);
    }

    centroids = malloc(sizeof(float *) * K);
    assert(centroids != NULL);

    for (i = 0; i < K; ++i) {
        centroids[i] = malloc(sizeof(float) * d);
        assert(centroids[i] != NULL);
    }

    sums = malloc(sizeof(float *) * K);
    assert(sums != NULL);

    for (i = 0; i < K; ++i) {
        sums[i] = malloc(sizeof(float) * d);
        assert(sums[i] != NULL);
    }

    clusters = malloc(N * sizeof(int));
    assert(clusters != NULL);

    csums = sums, cobs = obs, ccentroids = centroids;
    cclusters = clusters;
    algorithm(MAX_ITER, cobs, ccentroids, cclusters, csums);


    /* free memory */
    for (i = N - 1; i >= 0; --i) {
        free(obs[i]);
    }
    free(obs);
    for (i = K - 1; i >= 0; --i) {
        free(centroids[i]);
    }
    free(centroids);
    for (i = K - 1; i >= 0; --i) {
        free(sums[i]);
    }
    free(sums);
    free(clusters);
    return 0;

}

void initObs(float **obs, float **centroids) {
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


double norm(float *x, float *cluster) {
    double sum;
    int i;
    sum = 0;
    for (i = 0; i < d; ++i) {
        sum += (x[i] - cluster[i]) * (x[i] - cluster[i]);
    }
    return sum;
}

int assignCluster(float *x, float **centroids) {
    double sum;
    double tempSum;
    int minCluster;
    int i;
    sum = norm(x, centroids[0]);
    minCluster = 0;

    for (i = 1; i < K; ++i) {
        tempSum = norm(x, centroids[i]);
        if (tempSum < sum) {
            sum = tempSum;
            minCluster = i;
        }
    }
    return minCluster;
}

void assignAllObservations(float **obs, float **centroids, int *clusters) {
    int i;
    for (i = 0; i < N; ++i) {
        clusters[i] = assignCluster(obs[i], centroids);
    }
}

void resetSums(float **sums) {
    int i, j;
    for (i = 0; i < K; ++i)
        for (j = 0; j < d; ++j)
            sums[i][j] = 0;
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
char updateCentroids(float **obs, float **centroids, int *clusters, float **sums) {
    int *clusterSizes;
    char changedAny;
    float *tempCentroid;
    int i, j;
    resetSums(sums);
    changedAny = 0;
    clusterSizes = calloc(K, sizeof(int));
    tempCentroid = malloc(d * sizeof(float));
    for (i = 0; i < N; ++i) {
        for (j = 0; j < d; ++j) {
            sums[clusters[i]][j] += obs[i][j];
        }
        clusterSizes[clusters[i]]++;
    }
    for (i = 0; i < K; ++i) {
        for (j = 0; j < d; ++j) {
            if (clusterSizes[i] != 0)
                tempCentroid[j] = sums[i][j] / clusterSizes[i];
            else
                tempCentroid[j] = centroids[i][j];
            if (!changedAny && tempCentroid[j] != centroids[i][j])
                changedAny = 1;
            centroids[i][j] = tempCentroid[j];
        }
    }
    free(clusterSizes);
    free(tempCentroid);
    return changedAny;
}


void algorithm(int MAX_ITER, float **obs, float **centroids, int *clusters, float **sums) {
    char changedCluster;
    int i, j;
    changedCluster = 1;
    initObs(obs, centroids);

    for (i = 0; i < MAX_ITER && changedCluster; ++i) {
        assignAllObservations(obs, centroids, clusters);
        changedCluster = updateCentroids(obs, centroids, clusters, sums);
    }
    for (i = 0; i < K; ++i) {
        for (j = 0; j < d - 1; ++j) {
            printf("%.2f%c", centroids[i][j], ',');
        }
        printf("%.2f%c", centroids[i][d - 1], '\n');
    }

}
