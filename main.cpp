#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "tcm/tcm.hpp"
#include "channels/awgnchannel.hpp"
#include "modulators/modulator.hpp"

using namespace std;

vector<double> modeling(vector<complex<double>> &c, const vector<double> &snr_dB_arr, size_t num) {
    vector<double> Pe;
    AWGN_channel awgn;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(0, static_cast<int>(c.size()) - 1);

    for (auto &snr: snr_dB_arr) {
        cout << "SNR: " << snr << endl;
        auto sum = 0.0;
        size_t error_count = 0;
        for (auto n = 0; n < num; ++n) {
            auto i = distrib(gen);
            auto r = awgn.add_noise(c[i], snr);
            auto dmin = numeric_limits<double>::infinity();
            auto ind = 0;
            for (auto k = 0; k < c.size(); ++k) {
                double d = sqrt(pow(c[k].real() - r.real(), 2) + pow(c[k].imag() - r.imag(), 2));
                if (d < dmin) {
                    dmin = d;
                    ind = k;
                }
            }
            if (ind != i) {
                error_count++;
            }
        }
        Pe.push_back(static_cast<double>(error_count) / num);
    }
    return Pe;
}

int main() {
    //Time Test
//    vector<double> v (100000000);
//    AWGN_channel c;
//    auto start = std::chrono::system_clock::now();
//    v = c.add_noise(v,0, 1);
//    auto end = std::chrono::system_clock::now();
//
//    std::chrono::duration<double> diff = end - start;
//    std::cout << "Time to fill and iterate a vector of " << std::setw(9)
//              << v.size() << ": " << diff.count() << " s\n";
//
//    auto sum = 0.0;
//    for(auto &&vi: v){
//        sum += vi*vi;
//    }
//    cout << sum/ v.size() << endl;
//     params
    int q = 256;
    double T = 1.0 / 2400.0;
    double dt = T / 100.0;
    int f = 1800;
    // vector<double> f = {20.0/T, 19.0/T, 18.0/T, 17.0/T};
    double E = 1;
    //signals
    ofstream out("in.txt");
    auto arr = Modulator::qam(q, T, dt, f, E);

    int signal = 0;
    for(auto&& s : arr){
        double sum = 0;
        for (auto&& el : s){
            sum+= el*el;
        }
        sum *=dt;
        cout << setprecision(3) << "energy signal[" << signal++ << "]:" << sum << endl;
    }

    for (auto j = 0; j < arr[0].size(); ++j) {
        out << dt * j << ' ';
        for (auto &i: arr) {
            out << i[j] << ' ';
        }
        out << endl;
    }
    out.close();
    //constellation
    ofstream out2("constellation");
    auto foo = Modulator::qam_signal; //хранит указатель на функцию qam_signal
    auto constell = Modulator::constellation(q, T, dt, f, E, foo);
    for (int i = 0; i < q; ++i) {
        out2 << constell[i].real() << ' ' << constell[i].imag() << endl;
    }
    //awgn channel
    vector<double> snr_arr(25);
    int i = 15;
    for (auto &snr: snr_arr) {
        snr = i++;
    }
    ofstream out3("Pe");
    auto res = modeling(constell, snr_arr, 10000);
    for (auto j = 0; j < snr_arr.size(); ++j) {
        out3 << snr_arr[j] << ' ' << res[j] << ' ' << Modulator::fer_theor(pow(10, snr_arr[j] / 10), q, "qam") << endl;
    }

    //spectrum

    ofstream out4("Spectrum");
    auto sp = Modulator::spectrum(q, T, 1, f, E, "qam");
    for (int j = 0; j < sp[0].size(); ++j) {
        out4 << j << ' ' << sp[0][j] << endl;
    }

    vector<double> res_spectres;
    for (int j = 0; j < sp[0].size(); ++j) {
        res_spectres.push_back(sp[144][j] + sp[82][j] + sp[60][j]);
    }

    ofstream out_spectres("Spectrum_all");
    for (int j = 0; j < res_spectres.size(); ++j) {
        out_spectres << j << ' ' << res_spectres[j] << endl;
    }

    //for parseval
    int spectres = 0;
    for(auto&& sp : arr){
        double sum = 0;
        for (auto&& el : sp){
            sum+= el*el;
        }
        sum *=dt;
        cout << setprecision(3) << "energy spectrum[" << spectres++ << "]:" << sum << endl;
    }

    auto res2 = Modulator::grey_mapping(q,  "qam");
    for (int j = 0; j < res2.size(); ++j) {
        cout << res2[j] << ' ';
    }
    TCM tcm({{3,5}, {1,2}, {1,7}});
    vector<bool> msg ={1,0,1,0,1,1,1,0,1,0,0,1,0,1,1,1,1,0,0,0,0,0};
    auto enc_msg = tcm.encode(msg);
    for(const auto &b: enc_msg){
        cout << b?1:0;
    }
    cout << endl;
    for(const auto &b: msg){
        cout << b?1:0;
    }
    cout << endl;
    auto dec_msg = tcm.vitdec(enc_msg);
    for(const auto &b: dec_msg){
        cout << b?1:0;
    }
    return 0;
}

