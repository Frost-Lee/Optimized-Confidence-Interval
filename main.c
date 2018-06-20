#include <stdio.h>
#include <math.h>

#define DIVISION_PER_UNIT 10000
#define CHI_SQUARE_N 4
#define STEP_SIZE 1e-3

const static double gammaLnValues[20] =
        {0, -0.120782237635245, 0, 0.284682870472919, 0.693147180559945,1.20097360234707, 1.79175946922806,
         2.45373657084244, 3.17805383034795, 3.95781396761872, 4.78749174278204, 5.66256205985714, 6.5792512120101,
         7.53436423675873, 8.52516136106541, 9.549267257301, 10.6046029027452, 11.6893334207973, 12.8018274800815};

typedef struct {
    double lowerBound;
    double upperBound;
} Interval;

Interval getOptimizedInterval(double alpha);
double compositeSimpsonIntegration(double lowerBound, double upperBound);
double chiSquareDensityFunction(double x, double n);
double getGamma(double x);

int main() {
    Interval interval = getOptimizedInterval(0.05);
    printf("The interval is [%lf, %lf]", interval.lowerBound, interval.upperBound);
}

Interval getOptimizedInterval(double alpha) {
    double separationPoint = sqrt((CHI_SQUARE_N - 2) * M_E);
    double integration = compositeSimpsonIntegration(separationPoint - STEP_SIZE, STEP_SIZE);
    double frontPoint = separationPoint + STEP_SIZE;
    double rarePoint = separationPoint - STEP_SIZE;
    while (integration < 1.0 - alpha) {
        if (chiSquareDensityFunction(frontPoint, CHI_SQUARE_N) >
            chiSquareDensityFunction(rarePoint, CHI_SQUARE_N)) {
            integration = compositeSimpsonIntegration(rarePoint, frontPoint + STEP_SIZE);
            frontPoint += STEP_SIZE;
        } else {
            integration = compositeSimpsonIntegration(rarePoint, frontPoint + STEP_SIZE);
            rarePoint -= STEP_SIZE;
        }
    }
    Interval interval;
    interval.lowerBound = rarePoint;
    interval.upperBound = frontPoint;
    return interval;
}

double compositeSimpsonIntegration(double lowerBound, double upperBound) {
    double division = DIVISION_PER_UNIT * (upperBound - lowerBound);
    double steps = (upperBound - lowerBound) / division / 2.0;
    double integration = 0;
    double sumOfEven = 0;
    double sumOfOdd = 0;
    for (int i = 2; i < 2 * division; i += 2) {
        sumOfEven += chiSquareDensityFunction(lowerBound + i * steps, CHI_SQUARE_N);
    }
    for (int i = 1; i < 2 * division; i += 2) {
        sumOfOdd += chiSquareDensityFunction(lowerBound + i * steps, CHI_SQUARE_N);
    }
    integration = steps / 3.0 * (chiSquareDensityFunction(lowerBound, CHI_SQUARE_N) +
            4 * sumOfOdd + 2 * sumOfEven + chiSquareDensityFunction(upperBound, CHI_SQUARE_N));
    return integration;
}



// n is not expected to choose randomly
// Only values between 1 to 10, step 0.5 could be chosen from
double chiSquareDensityFunction(double x, double n) {
    if (x <= 0) {
        return 0;
    }
    double fractionPart = 1.0 / (pow(2, n / 2.0) * getGamma(n / 2.0));
    double plainPart = pow(x, n / 2.0 - 1) * pow(M_E, - x / 2.0);
    return fractionPart * plainPart;
}

// Not every value can be got
double getGamma(double x) {
    double gammaLn = gammaLnValues[(int)x * 2 - 2];
    double gammaValue = pow(M_E, gammaLn);
    return gammaValue;
}

