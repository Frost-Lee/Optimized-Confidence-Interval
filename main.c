#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DIVISION_PER_UNIT 10
#define INITIAL_STEP 0.1
#define STEP_SIZE 0.01
#define MAX_GAMMA 34

int CHI_SQUARE_N = 15;
int F_M = 21;
int F_N = 21;

const static double gammaLnValues[MAX_GAMMA + 1] =
        {0, -0.120782237635245, 0, 0.284682870472919, 0.693147180559945,1.20097360234707, 1.79175946922806,
         2.45373657084244, 3.17805383034795, 3.95781396761872, 4.78749174278204, 5.66256205985714, 6.5792512120101,
         7.53436423675873, 8.52516136106541, 9.549267257301, 10.6046029027452, 11.6893334207973, 12.8018274800815,
         13.9406252194038, 15.1044125730755, 16.2920004765672, 17.5023078458739, 18.7343475119364, 19.9872144956619,
         21.2600761562447, 22.5521638531234, 23.8627658416891, 25.1912211827387, 26.5369144911156, 27.8992713838409,
         29.2777545150408, 30.6718601060807};

typedef struct {
    double lowerBound;
    double upperBound;
} Interval;

typedef enum {
    false, true
} bool;

// Calculation related functions
Interval getOptimizedInterval(double alpha, double(*function)(double), double(*separation)());
double compositeSimpsonIntegration(double(*function)(double), double lowerBound, double upperBound);
double chiSquareDistribution(double x);
double fDistribution(double x);
double getChiSquareSeparation();
double getFSeparation();
double chiSquareDensityFunction(double x, double n);
double fDensityFunction(double x, double m, double n);
double getGamma(double x);

// Control related functions
char getDistributionChoice();
void getFDistributionParameter();
void getChiSquareDistributionParameter();
double getAlpha();



int main() {
//    double alphaValues[] = {0.02, 0.05, 0.1, 0.2};
//    for (int i = 0; i < 4; i ++) {
//        Interval interval = getOptimizedInterval(alphaValues[i], fDistribution, getFSeparation);
//        printf("[%.2lf, %.2lf]\t", interval.lowerBound, interval.upperBound);
//    }
    printf("This program can help you to get better optimized interval in F Distribution or"
           " Chi-Square Distribution.\n");
    while (true) {
        printf("Which distribution would you like to choose?\n");
        printf("'C' for Chi-Square distribution, 'F' for F distribution, and 'Q' for quit.\n");
        Interval interval;
        char choice = getDistributionChoice();
        double alpha = getAlpha();
        switch (choice) {
            case 'F':
                getFDistributionParameter();
                interval = getOptimizedInterval(alpha, fDistribution, getFSeparation);
                break;
            case 'C':
                getChiSquareDistributionParameter();
                interval = getOptimizedInterval(alpha, chiSquareDistribution, getChiSquareSeparation);
                break;
            case 'Q':
                exit(0);
            default:
                printf("This case shall not happen.\n");
                break;
        }
        printf("The interval is [%.2lf, %.2lf].\n\n", interval.lowerBound, interval.upperBound);
        getchar();
    }
}

Interval getOptimizedInterval(double alpha, double(*function)(double), double(*separation)()) {
    bool isBoundaryMode = false;
    double separationPoint = separation();
    double integration = compositeSimpsonIntegration(function, separationPoint - INITIAL_STEP / 2,
            separationPoint + INITIAL_STEP / 2);
    double frontPoint = separationPoint + INITIAL_STEP / 2;
    double rarePoint = separationPoint - INITIAL_STEP / 2;
    while (integration < 1.0 - alpha) {
        if (rarePoint <= 0) {
            isBoundaryMode = true;
        }
        if (isBoundaryMode) {
            rarePoint = 0;
            integration = compositeSimpsonIntegration(function, rarePoint, frontPoint + STEP_SIZE);
            frontPoint += STEP_SIZE;
        } else {
            double rareValue = function(rarePoint);
            double frontValue = function(frontPoint);
            if (frontValue > rareValue) {
                integration = compositeSimpsonIntegration(function, rarePoint, frontPoint + STEP_SIZE);
                frontPoint += STEP_SIZE;
            } else {
                integration = compositeSimpsonIntegration(function, rarePoint - STEP_SIZE, frontPoint);
                rarePoint -= STEP_SIZE;
            }
        }
    }
    Interval interval;
    interval.lowerBound = rarePoint;
    interval.upperBound = frontPoint;
    return interval;
}

double compositeSimpsonIntegration(double(*function)(double), double lowerBound, double upperBound) {
    double division = DIVISION_PER_UNIT * (upperBound - lowerBound);
    double steps = (upperBound - lowerBound) / division / 2.0;
    double integration = 0;
    double sumOfEven = 0;
    double sumOfOdd = 0;
    for (int i = 2; i < 2 * division; i += 2) {
        sumOfEven += function(lowerBound + i * steps);
    }
    for (int i = 1; i < 2 * division; i += 2) {
        sumOfOdd += function(lowerBound + i * steps);
    }
    integration = steps / 3.0 * (function(lowerBound) + 4 * sumOfOdd + 2 * sumOfEven + function(upperBound));
    return integration;
}

double chiSquareDistribution(double x) {
    return chiSquareDensityFunction(x, CHI_SQUARE_N);
}

double fDistribution(double x) {
    return fDensityFunction(x, F_M, F_N);
}

double getChiSquareSeparation() {
    return CHI_SQUARE_N - 2;
}

double getFSeparation() {
    return (F_M * F_N - 2.0 * F_N) / (F_M * F_N + 2.0 * F_M);
}

// The density function of Chi-Square distribution
// n is NOT expected to choose randomly
// Only values between 1 to 10, step 0.5 could be chosen from
double chiSquareDensityFunction(double x, double n) {
    if (x <= 0) {
        return 0;
    }
    const double part_1 = pow(x / 2.0, n / 2.0);
    const double part_2 = pow(M_E, - x / 2.0);
    const double part_3 = getGamma(n / 2.0);
    return part_1 * part_2 / x / part_3;
}

// The density function of F distribution
// n is NOT expected to choose randomly
// Only values between 1 to 10, step 0.5 could be chosen from
double fDensityFunction(double x, double m, double n) {
    if (x <= 0) {
        return 0;
    }
    const double part_1 = getGamma((m + n) / 2.0) / (getGamma(n / 2.0) * getGamma(m / 2.0));
    const double part_2 = pow(m, m / 2.0) * pow(n, n / 2.0) * pow(x, m / 2.0 - 1);
    const double part_3 = pow(n + m * x, -(m + n) / 2.0);
    return part_1 * part_2 * part_3;
}

// Not every value can be fetched
double getGamma(double x) {
    int slice = (int)(x * 2) - 2;
    if (slice > MAX_GAMMA || slice < 0) {
        printf("The gamma function value required in this question is not in the gamma function list.\n");
        printf("Try to add values to 'gammaLnValues' manually.\n");
        exit(1);
    }
    double gammaLn = gammaLnValues[slice];
    double gammaValue = pow(M_E, gammaLn);
    return gammaValue;
}

char getDistributionChoice() {
    char choice;
    do {
        scanf("%c", &choice);
        if (choice == 'C') {
            return 'C';
        } else if (choice == 'F') {
            return 'F';
        } else if (choice == 'Q') {
            printf("Thanks for using.\n");
            exit(0);
        } else {
            getchar();
            printf("Invalidate choice, try again.\n");
        }
    } while (true);
}

void getFDistributionParameter() {
    printf("Enter the two parameters (M and N, the order obeys the order of Beta function in F distribution)"
           " of F Distribution.\n");
    scanf("%d %d", &F_M, &F_N);
}

void getChiSquareDistributionParameter() {
    printf("Enter the parameter of Chi-Square Distribution.\n");
    scanf("%d", &CHI_SQUARE_N);
}

double getAlpha() {
    double alpha;
    printf("Now enter the confidence value (1 - alpha).\n");
    scanf("%lf", &alpha);
    return (1 - alpha) > 0 ? (1 - alpha) : 0.05;
}