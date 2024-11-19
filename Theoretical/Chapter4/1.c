#include <math.h>
#include <stdio.h>

int main() {
    float value = 1 - cosf(1.0 / 4.0);
    printf("cos(1/4) in single precision float: %f\n", value);
    return 0;
}