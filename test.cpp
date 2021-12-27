#include <iostream>
#include <fstream>
#include <sstream>

#include <chrono>

#include <cstdlib>

#ifdef __APPLE__
#include "OpenCL/cl.hpp"
#else
#include <CL/cl.hpp>
#endif
#include "clmp/clmp.h"
#include "ckks/context.h"
#include "ckks/secretkey.h"
#include "ckks/param.h"
#include "ckks/scheme.h"
#include "ckks/plaintext.h"
#include "ckks/ciphertext.h"
#include "ckks/key.h"
/*
std::string read_kernel_code_from_file(std::string _path) {
    std::stringstream buf;
    std::ifstream fin(_path);

    buf << fin.rdbuf();

    return buf.str();
}
*/
int main() {
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    if (platforms.size()==0) {
        std::cout<<" No platforms found. Check OpenCL installation!" << std::endl;;
        exit(1);
    }

    std::cout << platforms.size() << " platforms are found!." << std::endl;

    // Detect OpenCL Device
    cl::Platform platform = platforms[0];
    std::cout << "Using platform: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;

    std::vector<cl::Device> devices;
    platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);

    if(devices.size() == 0){
        std::cout<<" No devices found. Check OpenCL installation!" << std::endl;;
        exit(1);
    }

    std::cout << devices.size() << " devices are found!." << std::endl;

    cl::Device device=devices[2];
    std::cout<< "Using device: "<<device.getInfo<CL_DEVICE_NAME>() << std::endl;

    cl::Context cl_ctx({device});


    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    CKKS::context ckks_ctx(&cl_ctx);

    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::cout << "boot time : " << sec.count() << " seconds" << std::endl;


    CKKS::secretkey s_key(ckks_ctx);
    CKKS::scheme ckks_scheme(ckks_ctx, s_key);

    clmp_uint64 log_n = 3;
    clmp_uint64 n = 1 << log_n;
    clmp_uint64 log_p = 30;
    clmp_uint64 log_q = 100;
    clmp_uint64 num = 1;

    std::complex<double> *input_1 = new std::complex<double>[n * num];

    for (clmp_uint64 i = 0; i < n * num; i++) {
        input_1[i] = i;
    }

    for (clmp_uint64 i = 0; i < n * num; i++) {
        std::cout << input_1[i] << ", ";
    }
    std::cout << std::endl << std::endl;

    CKKS::ciphertext cipher_1;
    start = std::chrono::system_clock::now();
    ckks_scheme.encrypt(cipher_1, input_1, n, num, log_p, log_q);
    sec = std::chrono::system_clock::now() - start;
    std::cout << "encryption time : " << sec.count() << " seconds" << std::endl;


    start = std::chrono::system_clock::now();
    std::complex<double> *output_1 = ckks_scheme.decrypt(s_key, cipher_1);
    sec = std::chrono::system_clock::now() - start;
    std::cout << "decryption time : " << sec.count() << " seconds" << std::endl;

    for (clmp_uint64 i = 0; i < n * num; i++) {
        std::cout << output_1[i] << ", ";
    }
    std::cout << std::endl << std::endl;

    delete [] input_1;
    delete [] output_1;

    return 0;
}
