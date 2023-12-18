#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <windows.h>

typedef struct {
    double real;
    double imag;
} complex;


typedef struct {
    double resistance[2];
    double capacitance;
    double inductance;
    double currentFrequency;
    double angularFrequency;
} RLC_Circuit;


void calculateRLC1(complex *Numerator, complex *Denominator, RLC_Circuit circuit1);

void text();

void calculateRLC2(complex *Numerator, complex *Denominator, RLC_Circuit circuit2);

void calculateRLC3(complex *Numerator, complex *Denominator, RLC_Circuit circuit3);

void calculateRLC4(complex *Numerator, complex *Denominator, RLC_Circuit circuit4);

void circuitParameterInput(RLC_Circuit *, int, int (*)(double));

void calculateInRange(RLC_Circuit *circuit, void (*)(complex *, complex *, RLC_Circuit));

void choosingCircuit(int *);

double input(int(*)(double));

void rangeInput(double *, double *, double *);

int isValidParameter(double);

complex calculateReactiveResistance1(void (*f)(complex *, complex *, RLC_Circuit), RLC_Circuit circuit);

complex calculateRatio(complex, complex);


int main() {
    system("chcp 65001");
    system("cls");
    text();
    char a[20];
    do {
        RLC_Circuit newCircuit;
        int chosenCircuit = 0;
        void (*c[4])(complex *, complex *, RLC_Circuit) = {calculateRLC1, calculateRLC2, calculateRLC3, calculateRLC4};
        choosingCircuit(&chosenCircuit);
        circuitParameterInput(&newCircuit, chosenCircuit, isValidParameter);
        calculateInRange(&newCircuit, c[chosenCircuit - 1]);
    } while (getch() != 27);
}

complex calculateReactiveResistance1(void (*f)(complex *, complex *, RLC_Circuit), RLC_Circuit circuit) {
    complex Numerator, Denominator;
    f(&Numerator, &Denominator, circuit);
    return calculateRatio(Numerator, Denominator);
}


void calculateRLC1(complex *Numerator, complex *Denominator, RLC_Circuit circuit1) {
    Numerator->real = circuit1.inductance / circuit1.capacitance;
    Numerator->imag = -circuit1.resistance[0] / (circuit1.angularFrequency * circuit1.capacitance);
    Denominator->real = circuit1.resistance[0];
    Denominator->imag = (circuit1.angularFrequency * circuit1.inductance -
                         (1.0 / (circuit1.angularFrequency * circuit1.capacitance)));
}

void calculateRLC2(complex *Numerator, complex *Denominator, RLC_Circuit circuit2) {
    Numerator->real = circuit2.inductance / circuit2.capacitance;
    Numerator->imag = circuit2.resistance[0] / (circuit2.angularFrequency * circuit2.capacitance);
    Denominator->real = circuit2.resistance[0];
    Denominator->imag = (circuit2.angularFrequency * circuit2.inductance -
                         (1.0 / (circuit2.angularFrequency * circuit2.capacitance)));
}

void calculateRLC3(complex *Numerator, complex *Denominator, RLC_Circuit circuit3) {
    Numerator->real = circuit3.resistance[0] * circuit3.resistance[1];
    Numerator->imag = circuit3.resistance[0] *
                      (circuit3.angularFrequency * circuit3.inductance -
                       (1.0 / (circuit3.angularFrequency * circuit3.capacitance)));
    Denominator->real = circuit3.resistance[0] + circuit3.resistance[1];
    Denominator->imag = (circuit3.angularFrequency * circuit3.inductance -
                         (1.0 / (circuit3.angularFrequency * circuit3.capacitance)));
}

void calculateRLC4(complex *Numerator, complex *Denominator, RLC_Circuit circuit4) {
    Numerator->real = circuit4.resistance[0] * circuit4.resistance[1] + (circuit4.inductance / circuit4.capacitance);
    Numerator->imag = circuit4.angularFrequency * circuit4.inductance * circuit4.resistance[0] -
                      (circuit4.resistance[1] / (circuit4.angularFrequency * circuit4.capacitance));
    Denominator->real = circuit4.resistance[0] + circuit4.resistance[1];
    Denominator->imag = (circuit4.angularFrequency * circuit4.inductance -
                         (1.0 / (circuit4.angularFrequency * circuit4.capacitance)));
}

complex calculateRatio(complex num, complex den) {
    complex Z;
    Z.real = (num.real * den.real + num.imag * den.imag) / (den.real * den.real + den.imag * den.imag);
    Z.imag = (num.imag * den.real - num.real * den.imag) / (den.real * den.real + den.imag * den.imag);
    return Z;
}

void calculateInRange(RLC_Circuit *circuit, void (*f)(complex *, complex *, RLC_Circuit)) {
    double startFreq = 0, df = 0, finishFrequency = 0;
    rangeInput(&startFreq, &finishFrequency, &df);
    int len = (int) (fabs(startFreq - finishFrequency) / df) + 1;
    complex Z;
    circuit->currentFrequency = startFreq;
    circuit->angularFrequency = 2.0 * M_PI * startFreq;
    for (int i = 0; i < len; i++) {
        Z = calculateReactiveResistance1(f, *circuit);
        printf("freq: %lf\t", circuit->currentFrequency);
        printf("%.15lf %+.15lfi\t", Z.real, Z.imag);
        printf("%.15lf\n", 1.0 / sqrt(circuit->inductance * circuit->capacitance));
        circuit->currentFrequency += df;
        circuit->angularFrequency += 2.0 * M_PI * df;
    }
}

void choosingCircuit(int *a) {
    printf("Choose a circuit");
    while (scanf("%d", a) != 1 || getchar() != '\n' || *a > 4 || *a < 1) {
        fflush(stdin);
        printf("Invalid input! Enter a positive int number!");
    }
}

void circuitParameterInput(RLC_Circuit *circuit, int circuitNumber, int (*isValidParameter)(double)) {
    printf("Enter R1:");
    circuit->resistance[0] = input(isValidParameter);
    if (circuitNumber != 1 && circuitNumber != 2) {
        printf("Enter R2:");
        circuit->resistance[1] = input(isValidParameter);
    }
    printf("Enter C:");
    circuit->capacitance = input(isValidParameter);
    printf("Enter L:");
    circuit->inductance = input(isValidParameter);
}

double input(int (*f)(double a)) {
    double a = 0;
    while (scanf("%lf", &a) != 1 || getchar() != '\n' || f(a) == 0) {
        fflush(stdin);
        printf("Invalid input! Enter a positive real number!");
    }
    return a;
}

int isValidParameter(double a) {
    if (a < 0 || a > 1000) return 0;
    return 1;
}


void rangeInput(double *f1, double *f2, double *df) {
    printf("Enter start frequency:");
    *f1 = input(isValidParameter);
    printf("Enter finish frequency:");
    *f2 = input(isValidParameter);
    printf("Enter step:");
    *df = input(isValidParameter);
    if (*f1 > *f2) *df *= -1.0;

}

void text() {
    printf("This program calculates the reactive resistance of one of 4 RLC circuits"
           "The first circuit contains 1 resistor, 1 inductor and 1 capacitor"
           "The second circuit is like the first one but capacitor and inductor are swapped"
           "The third one and the fourth both have 2 resistor. Take a look at the circuits: ");

    printf("\nThe first:\n");
    printf("●—┳——███—————∩∩∩∩—┳—●\n"
           "  ┃   R       L   ┃\n"
           "  ┃               ┃\n"
           "  ┃        C      ┃\n"
           "  ┖————————││—————┚\n\n");

    printf("The second:\n");
    printf("\n●—┳——███—————││———┳—●\n"
           "  ┃   R      C    ┃\n"
           "  ┃               ┃\n"
           "  ┃        L      ┃\n"
           "  ┖———————∩∩∩∩————┚\n\n");
    printf("The third:\n");
    printf("●—┳——███—————││———┓\n"
           "  ┃  R2      C    ┃\n"
           "  ┃               ┃\n"
           "  █ R1     L      ┃\n"
           "●—┸———————∩∩∩∩————┚\n\n");

    printf("The fourth:\n");
    printf("●—┳—————███———————┓\n"
           "  █ R1   R2       ┃\n"
           "  ┻               ┃\n"
           "  ┳        L      ┃\n"
           "●—┸———————∩∩∩∩————┚\n\n");
}
