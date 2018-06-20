#include <stdio.h>
#include <math.h>

#define DIVISION_PER_UNIT 100
#define CHI_SQUARE_N 30
#define INITIAL_STEP 0.5
#define STEP_SIZE 0.01

const static double gammaLnValues[35] =
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

Interval getOptimizedInterval(double alpha);
double compositeSimpsonIntegration(double lowerBound, double upperBound);
double chiSquareDensityFunction(double x, double n);
double getGamma(double x);

int main() {
    Interval interval = getOptimizedInterval(0.05);
    printf("The interval is [%lf, %lf].\n", interval.lowerBound, interval.upperBound);
}

Interval getOptimizedInterval(double alpha) {
    double separationPoint = CHI_SQUARE_N - 2;
    double integration = compositeSimpsonIntegration(separationPoint - INITIAL_STEP / 2,
            separationPoint + INITIAL_STEP / 2);
    double frontPoint = separationPoint + INITIAL_STEP / 2;
    double rarePoint = separationPoint - INITIAL_STEP / 2;
    while (integration < 1.0 - alpha) {
        double rareValue = chiSquareDensityFunction(rarePoint, CHI_SQUARE_N);
        double frontValue = chiSquareDensityFunction(frontPoint, CHI_SQUARE_N);
        if (frontValue > rareValue) {
            integration = compositeSimpsonIntegration(rarePoint, frontPoint + STEP_SIZE);
            frontPoint += STEP_SIZE;
        } else {
            integration = compositeSimpsonIntegration(rarePoint - STEP_SIZE, frontPoint);
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
    const double part_1 = pow(x / 2.0, n / 2.0);
    const double part_2 = pow(M_E, - x / 2.0);
    const double part_3 = getGamma(n / 2.0);
    return part_1 * part_2 / x / part_3;
}

// Not every value can be fetched
double getGamma(double x) {
    double gammaLn = gammaLnValues[(int)x * 2 - 2];
    double gammaValue = pow(M_E, gammaLn);
    return gammaValue;
}

