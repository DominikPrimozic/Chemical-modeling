#include <FFT.h>

const double PI = acos(-1);
void dfft(std::vector<std::complex<double>>& a){
    int n = a.size();
    if (n == 1)
        return;
    std::vector<std::complex<double>> a0(n / 2), a1(n / 2);
    for (int i = 0; 2 * i < n; i++) {
        a0[i] = a[2*i];
        a1[i] = a[2*i+1];
    }
    dfft(a0);
    dfft(a1);

    double ang = 2 * PI / n;
    std::complex<double> w(1), wn(cos(ang), sin(ang));
    for (int i = 0; 2 * i < n; i++) {
        a[i] = a0[i] + w * a1[i];
        a[i + n / 2] = a0[i] - w * a1[i];
        w *= wn;
    }
}