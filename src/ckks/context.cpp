#include "ckks/context.h"

#include <complex>

#include <chrono>

std::string read_kernel_code_from_file(std::string _path) {
    std::stringstream buf;
    std::ifstream fin(_path);

    buf << fin.rdbuf();

    return buf.str();
}

namespace CKKS {
    context::context(const cl::Context *_cl_ctx) {
        // OpenCL context check...
        if (_cl_ctx == NULL) {
            std::cout << "OpenCL context is not specified...." << std::endl;
            std::cout << "Use default context..." << std::endl;

            std::vector<cl::Platform> platforms;
            cl::Platform::get(&platforms);

            if (platforms.size()==0) {
                std::cout<<" No platforms found. Check OpenCL installation!" << std::endl;;
                exit(1);
            }

            // Detect OpenCL Device
            cl::Platform platform = platforms[0];
            std::cout << "Using platform: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;

            std::vector<cl::Device> devices;
            platform.getDevices(CL_DEVICE_TYPE_CPU | CL_DEVICE_TYPE_GPU, &devices);

            if(devices.size()==0){
                std::cout<<" No devices found. Check OpenCL installation!" << std::endl;;
                exit(1);
            }

            cl::Device device=devices[0];
            std::cout<< "Using device: "<<device.getInfo<CL_DEVICE_NAME>() << std::endl;

            cl_ctx = new cl::Context({device});
        }
        else {
            cl_ctx = new cl::Context(*_cl_ctx);
        }

        cl::Program::Sources sources;
        std::string clmp_src = read_kernel_code_from_file("./src/clmp.cl");
        sources.push_back({clmp_src.c_str(), clmp_src.length()});

        std::string ckks_src = read_kernel_code_from_file("./src/ckks.cl");
        sources.push_back({ckks_src.c_str(), ckks_src.length()});

        cl_program = cl::Program(*cl_ctx, sources);

        if(cl_program.build({cl_ctx->getInfo<CL_CONTEXT_DEVICES>()[0]}, "-I ./include")!=CL_SUCCESS){
            std::cout << " Error building: " << cl_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(cl_ctx->getInfo<CL_CONTEXT_DEVICES>()[0]) << std::endl;
            exit(1);
        }

        cl_queue = cl::CommandQueue(*cl_ctx, cl_ctx->getInfo<CL_CONTEXT_DEVICES>()[0]);

        clmpz_set_ui(Q, 1);
        clmpz_set_ui(QQ, 1);
        clmpz_shift_left_ui(Q, Q, logQ);
        clmpz_shift_left_ui(QQ, QQ, logQQ);

        //std::cout << clmpz_to_string(QQ) << std::endl;

        // Ring Generation...
        cl_q_pows = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * (logQQ + 1));

        clmpz_t tmp;
        clmpz_set_i(tmp, 1);

        clmpz_set_ui(q_pows[0], 1);

        q_pows_k_d[0] = clmpz_bits(q_pows[0]) * 2;

        clmpz_shift_left_ui(q_pows_r[0], tmp, q_pows_k_d[0]);
        clmpz_div(q_pows_r[0], q_pows_r[0], q_pows[0]);
        for (clmp_uint64 i = 1; i < logQQ + 1; i++) {
            clmpz_shift_left_ui(q_pows[i], q_pows[i - 1], 1);

            q_pows_k_d[i] = clmpz_bits(q_pows[i]) * 2;

            clmpz_shift_left_ui(q_pows_r[i], tmp, q_pows_k_d[i]);
            clmpz_div(q_pows_r[i], q_pows_r[i], q_pows[i]);
        }

        cl_queue.enqueueWriteBuffer(cl_q_pows, CL_TRUE, 0, sizeof(clmpz_t) * (logQQ + 1), q_pows);


        cl_rot_group = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * Nh);

        clmp_uint64 five_pows = 1;
        for (clmp_uint64 i = 0; i < Nh; i++) {
            rot_group[i] = five_pows;
            five_pows *= 5;
            five_pows %= M;
        }

        cl_queue.enqueueWriteBuffer(cl_rot_group, CL_TRUE, 0, sizeof(long) * Nh, rot_group);


        cl_ksi_pows = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(std::complex<double>) * (M + 1));
        for (clmp_uint64 i = 0; i < M; i++) {
            double angle = 2.0 * M_PI * i / M;
            ksi_pows[i].real(cos(angle));
            ksi_pows[i].imag(sin(angle));
        }
        ksi_pows[M] = ksi_pows[0];

        cl_queue.enqueueWriteBuffer(cl_ksi_pows, CL_TRUE, 0, sizeof(std::complex<double>) * (M + 1), ksi_pows);


        cl_p_vec = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * nprimes);
        cl_p_r_vec = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * nprimes);
        cl_p_inv_vec = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * nprimes);
        cl_scaled_root_pows = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * nprimes * N);
        cl_scaled_root_inv_pows = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * nprimes * N);
        cl_scaled_N_inv = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * nprimes);

        clmp_uint64 prime = (1ULL << pbnd) + 1;
        for (clmp_uint64 i = 0; i < nprimes; i++) {
            while (true) {
                prime += M;

                if (is_prime(prime, 200)) {
                    p_vec[i] = prime;
                    break;
                }
            }
        }

        for (clmp_uint64 i = 0; i < nprimes; i++) {
            p_r_vec[i] = (static_cast<unsigned __int128>(1) << kbar2) / p_vec[i];
            p_inv_vec[i] = inv(p_vec[i]);

            clmp_uint64 root = find_Nth_root_of_unity(M, p_vec[i]);
            clmp_uint64 root_inv = inv_mod(root, p_vec[i]);

            clmp_uint64 N_inv = inv_mod(N, p_vec[i]);
            scaled_N_inv[i] = mul_mod(N_inv, (1ULL << 32), p_vec[i]);
            scaled_N_inv[i] = mul_mod(scaled_N_inv[i], (1ULL << 32), p_vec[i]);

            clmp_uint64 power = 1;
            clmp_uint64 power_inv = 1;

            for (clmp_uint64 j = 0; j < N; j++) {
                clmp_uint32 j_prime = reverse_bit(static_cast<clmp_uint32>(j)) >> (32 - logN);

                clmp_uint64 root_pow = power;
                scaled_root_pows[i * N + j_prime] = mul_mod(root_pow, (1ULL << 32), p_vec[i]);
                scaled_root_pows[i * N + j_prime] = mul_mod(scaled_root_pows[i * N + j_prime], (1ULL << 32), p_vec[i]);
                std::cout << scaled_root_pows[i * N + j_prime] << std::endl;

                clmp_uint64 root_pow_inv = power_inv;
                scaled_root_inv_pows[i * N + j_prime] = mul_mod(root_pow_inv, (1ULL << 32), p_vec[i]);
                scaled_root_inv_pows[i * N + j_prime] = mul_mod(scaled_root_inv_pows[i * N + j_prime], (1ULL << 32), p_vec[i]);

                power = mul_mod(power, root, p_vec[i]);
                power_inv = mul_mod(power_inv, root_inv, p_vec[i]);
            }
        }

        cl_queue.enqueueWriteBuffer(cl_p_vec, CL_TRUE, 0, sizeof(clmp_uint64) * nprimes, p_vec);
        cl_queue.enqueueWriteBuffer(cl_p_r_vec, CL_TRUE, 0, sizeof(clmp_uint64) * nprimes, p_r_vec);
        cl_queue.enqueueWriteBuffer(cl_p_inv_vec, CL_TRUE, 0, sizeof(clmp_uint64) * nprimes, p_inv_vec);
        cl_queue.enqueueWriteBuffer(cl_scaled_root_pows, CL_TRUE, 0, sizeof(clmp_uint64) * nprimes * N, scaled_root_pows);
        cl_queue.enqueueWriteBuffer(cl_scaled_root_inv_pows, CL_TRUE, 0, sizeof(clmp_uint64) * nprimes * N, scaled_root_inv_pows);
        cl_queue.enqueueWriteBuffer(cl_scaled_N_inv, CL_TRUE, 0, sizeof(clmp_uint64) * nprimes, scaled_N_inv);


        cl_p_prod = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * nprimes);
        cl_p_prod_h = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * nprimes);
        cl_p_hat = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * nprimes * nprimes);
        cl_p_hat_inv_mod_p = cl::Buffer(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * nprimes * nprimes);

        for (clmp_uint64 i = 0; i < nprimes; i++) {
            if (i == 0) {
                clmpz_set_ul(p_prod[i], p_vec[i]);
            }
            else {
                clmpz_mul_ul(p_prod[i], p_prod[i - 1], p_vec[i]);
            }

            clmpz_div_ui(p_prod_h[i], p_prod[i], 2);

            for (clmp_uint64 j = 0; j < i + 1; j++) {
                clmpz_set_i(p_hat[i * nprimes + j], 1);

                for (clmp_uint64 k = 0; k < j; k++) {
                    clmpz_mul_ul(p_hat[i * nprimes + j], p_hat[i * nprimes + j], p_vec[k]);
                }

                for (clmp_uint64 k = j + 1; k < i + 1; k++) {
                    clmpz_mul_ul(p_hat[i * nprimes + j], p_hat[i * nprimes + j], p_vec[k]);
                }

                clmpz_t tmp;
                clmpz_mod_ul(tmp, p_hat[i * nprimes + j], p_vec[j]);
                p_hat_inv_mod_p[i * nprimes + j] = clmpz_to_ul(tmp);
                p_hat_inv_mod_p[i * nprimes + j] = inv_mod(p_hat_inv_mod_p[i * nprimes + j], p_vec[j]);


            }
        }

        cl_queue.enqueueWriteBuffer(cl_p_prod, CL_TRUE, 0, sizeof(clmpz_t) * nprimes, p_prod);
        cl_queue.enqueueWriteBuffer(cl_p_prod_h, CL_TRUE, 0, sizeof(clmpz_t) * nprimes, p_prod_h);
        cl_queue.enqueueWriteBuffer(cl_p_hat, CL_TRUE, 0, sizeof(clmpz_t) * nprimes * nprimes, p_hat);
        cl_queue.enqueueWriteBuffer(cl_p_hat_inv_mod_p, CL_TRUE, 0, sizeof(clmp_uint64) * nprimes * nprimes, p_hat_inv_mod_p);
    }

    context::~context() {
        delete cl_ctx;

        delete [] q_pows;
        delete [] q_pows_r;
        delete [] q_pows_k_d;
        delete [] rot_group;
        delete [] ksi_pows;

        delete [] p_vec;
        delete [] p_r_vec;
        delete [] p_inv_vec;
        delete [] scaled_root_pows;
        delete [] scaled_root_inv_pows;
        delete [] scaled_N_inv;

        delete [] p_prod;
        delete [] p_prod_h;
        delete [] p_hat;
        delete [] p_hat_inv_mod_p;
    }

    clmp_uint32 context::reverse_bit(const clmp_uint32 &_rhs) {
        clmp_uint32 res = _rhs;

        res = (((res & 0xaaaaaaaa) >> 1) | ((res & 0x55555555) << 1));
        res = (((res & 0xcccccccc) >> 2) | ((res & 0x33333333) << 2));
        res = (((res & 0xf0f0f0f0) >> 4) | ((res & 0x0f0f0f0f) << 4));
        res = (((res & 0xff00ff00) >> 8) | ((res & 0x00ff00ff) << 8));

        return ((res >> 16) | (res << 16));
    }

    clmp_uint64 context::pow(const clmp_uint64 &_lhs, const clmp_uint64 &_rhs) {
        clmp_uint64 res = 1;
        clmp_uint64 base = _lhs, exp = _rhs;

        while (exp > 0) {
            if (exp & 1) {
                res *= base;
            }
            exp >>= 1;
            base *= base;
        }

        return res;
    }

    clmp_uint64 context::inv(const clmp_uint64 &_rhs) {
        return pow(_rhs, static_cast<clmp_uint64>(-1));
    }

    clmp_uint64 context::mul_mod(const clmp_uint64 &_lhs, const clmp_uint64 &_rhs, const clmp_uint64 &_mod) {
        unsigned __int128 mul = static_cast<unsigned __int128>(_lhs) * _rhs;
    	mul %= static_cast<unsigned __int128>(_mod);
    	return static_cast<clmp_uint64>(mul);
    }

    clmp_uint64 context::pow_mod(const clmp_uint64 &_lhs, const clmp_uint64 &_rhs, const clmp_uint64 &_mod) {
        clmp_uint64 res = 1;
        clmp_uint64 base = _lhs, exp = _rhs;

        while (exp > 0)  {
            if (exp & 1) {
                res = mul_mod(res, base, _mod);
            }
            exp >>= 1;
            base = mul_mod(base, base, _mod);
        }
        return res;
    }

    clmp_uint64 context::inv_mod(const clmp_uint64 &_rhs, const clmp_uint64 &_mod) {
        return pow_mod(_rhs, _mod-2, _mod);
    }

    void context::find_prime_factors(std::vector<clmp_uint64> &_res, clmp_uint64 _rhs) {
        while (_rhs % 2  == 0) {
            _res.push_back(2);
            _rhs /= 2;
        }

        for (clmp_uint64 i = 3; i < sqrt(_rhs); i++) {
            while (_rhs % i == 0) {
                _res.push_back(i);
                _rhs /= i;
            }
        }

        if (_rhs > 2) {
            _res.push_back(_rhs);
        }
    }

    clmp_uint64 context::find_primitive_root(const clmp_uint64 &_mod) {
        std::vector<clmp_uint64> prime_factors;
        clmp_uint64 phi = _mod -1;

        find_prime_factors(prime_factors, phi);

        for (clmp_uint64 r = 2; r <= phi; r++) {
            bool flag = false;

            for (auto factor_iter = prime_factors.begin(); factor_iter != prime_factors.end(); factor_iter++) {
                if (pow_mod(r, phi/(*factor_iter), _mod) == 1) {
                    flag = true;
                    break;
                }
            }
            if (!flag) {
                return r;
            }
        }

        return 0;
    }

    clmp_uint64 context::find_Nth_root_of_unity(const long &_N, const clmp_uint64 &_mod) {
        clmp_uint64 res;

        res = find_primitive_root(_mod);

        if ((_mod - 1) % _N == 0) {
            clmp_uint64 factor = (_mod - 1) / _N;

            res = pow_mod(res, factor, _mod);

            return res;
        }
        else {
            return 0;
        }
    }

    bool context::is_prime(const clmp_uint64 &_rhs, const clmp_size_t &_n_iter) {
        if (_rhs < 2) {
            return false;
        }

        if ((_rhs != 2) && (_rhs % 2 == 0)) {
            return false;
        }

        clmp_uint64 s = _rhs - 1;
        while(s % 2 == 0) {
            s = s / 2;
        }

        for (size_t i = 0; i < _n_iter; i++) {
            clmp_uint64 rand_64 = static_cast<clmp_uint64>(rand()) << 32 | rand();
            clmp_uint64 temp = s;
            clmp_uint64 a = rand_64 % (_rhs - 1) + 1;
            clmp_uint64 mod = pow_mod(a, temp, _rhs);

            while((temp != _rhs - 1) && (mod != 1 ) && (mod != _rhs -1)) {
                mod = mul_mod(mod, mod, _rhs);
                temp = temp * 2;
            }
            if ((mod != _rhs - 1) && (temp % 2 == 0)) {
                return false;
            }
        }

        return true;
    }

    void context::fft_special(std::complex<double> *_vals, const long _len) {
        cl::Buffer cl_vals(*cl_ctx, CL_MEM_READ_WRITE, sizeof(std::complex<double>) * _len);

        cl_queue.enqueueWriteBuffer(cl_vals, CL_TRUE, 0, sizeof(std::complex<double>) * _len, _vals);

        cl::Kernel cl_fft_special(cl_program, "fft_special");
        cl_fft_special.setArg(0, cl_vals);
        cl_fft_special.setArg(1, _len);
        cl_fft_special.setArg(2, M);
        cl_fft_special.setArg(3, cl_rot_group);
        cl_fft_special.setArg(4, cl_ksi_pows);

        cl_queue.enqueueNDRangeKernel(cl_fft_special, cl::NullRange, cl::NDRange(1), cl::NullRange);

        cl_queue.enqueueReadBuffer(cl_vals, CL_TRUE, 0, sizeof(std::complex<double>) * _len, _vals);
    }

    void context::fft_special_inv(std::complex<double> *_vals, const long _len) {
        cl::Buffer cl_vals(*cl_ctx, CL_MEM_READ_WRITE, sizeof(std::complex<double>) * _len);

        cl_queue.enqueueWriteBuffer(cl_vals, CL_TRUE, 0, sizeof(std::complex<double>) * _len, _vals);

        cl::Kernel cl_fft_special_inv(cl_program, "fft_special_inv");
        cl_fft_special_inv.setArg(0, cl_vals);
        cl_fft_special_inv.setArg(1, _len);
        cl_fft_special_inv.setArg(2, M);
        cl_fft_special_inv.setArg(3, cl_rot_group);
        cl_fft_special_inv.setArg(4, cl_ksi_pows);

        cl_queue.enqueueNDRangeKernel(cl_fft_special_inv, cl::NullRange, cl::NDRange(1), cl::NullRange);

        cl_queue.enqueueReadBuffer(cl_vals, CL_TRUE, 0, sizeof(std::complex<double>) * _len, _vals);
    }

    void context::encode(clmpz_t *_res, const std::complex<double> *_vals, const long _len, const long _num, const long _log_p) {
        cl::Buffer cl_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _num);
        cl::Buffer cl_vals(*cl_ctx, CL_MEM_READ_WRITE, sizeof(std::complex<double>) * _len * _num);

        cl_queue.enqueueWriteBuffer(cl_vals, CL_TRUE, 0, sizeof(std::complex<double>) * _len * _num, _vals);

        cl::Kernel cl_encode(cl_program, "encode");
        cl_encode.setArg(0, cl_res);
        cl_encode.setArg(1, cl_vals);
        cl_encode.setArg(2, _len);
        cl_encode.setArg(3, _log_p);
        cl_encode.setArg(4, Nh);
        cl_encode.setArg(5, M);
        cl_encode.setArg(6, cl_rot_group);
        cl_encode.setArg(7, cl_ksi_pows);

        cl_queue.enqueueNDRangeKernel(cl_encode, cl::NullRange, cl::NDRange(_num), cl::NullRange);

        cl_queue.enqueueReadBuffer(cl_res, CL_TRUE, 0, sizeof(clmpz_t) * N * _num, _res);
    }

    void context::decode(std::complex<double> *_res, const clmpz_t *_vals, const long _len, const long _num, const long _log_p, const long _log_q) {
        cl::Buffer cl_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(std::complex<double>) * _len * _num);
        cl::Buffer cl_vals(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _num);

        cl_queue.enqueueWriteBuffer(cl_vals, CL_TRUE, 0, sizeof(clmpz_t) * N * _num, _vals);

        cl::Kernel cl_decode(cl_program, "decode");
        cl_decode.setArg(0, cl_res);
        cl_decode.setArg(1, cl_vals);
        cl_decode.setArg(2, _len);
        cl_decode.setArg(3, _log_p);
        cl_decode.setArg(4, _log_q);
        cl_decode.setArg(5, Nh);
        cl_decode.setArg(6, M);
        cl_decode.setArg(7, cl_rot_group);
        cl_decode.setArg(8, cl_ksi_pows);
        cl_decode.setArg(9, cl_q_pows);

        cl_queue.enqueueNDRangeKernel(cl_decode, cl::NullRange, cl::NDRange(_num), cl::NullRange);

        cl_queue.enqueueReadBuffer(cl_res, CL_TRUE, 0, sizeof(std::complex<double>) * _len * _num, _res);
    }

    void context::NTT(clmp_uint64 *_vals) {
        cl::Buffer cl_vals(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * N);

        cl_queue.enqueueWriteBuffer(cl_vals, CL_TRUE, 0, sizeof(clmp_uint64) * N, _vals);

        cl::Kernel cl_ntt(cl_program, "NTT");
        cl_ntt.setArg(0, cl_vals);
        cl_ntt.setArg(1, p_vec[0]);
        cl_ntt.setArg(2, p_inv_vec[0]);
        cl_ntt.setArg(3, cl_scaled_root_pows);
        cl_ntt.setArg(4, N);
        cl_ntt.setArg(5, logN);

        cl_queue.enqueueNDRangeKernel(cl_ntt, cl::NullRange, cl::NDRange(1), cl::NullRange);

        cl_queue.enqueueReadBuffer(cl_vals, CL_TRUE, 0, sizeof(clmp_uint64) * N, _vals);
    }

    void context::iNTT(clmp_uint64 *_vals) {
        cl::Buffer cl_vals(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * N);

        cl_queue.enqueueWriteBuffer(cl_vals, CL_TRUE, 0, sizeof(clmp_uint64) * N, _vals);

        cl::Kernel cl_intt(cl_program, "iNTT");
        cl_intt.setArg(0, cl_vals);
        cl_intt.setArg(1, p_vec[0]);
        cl_intt.setArg(2, p_inv_vec[0]);
        cl_intt.setArg(3, cl_scaled_root_inv_pows);
        cl_intt.setArg(4, scaled_N_inv[0]);
        cl_intt.setArg(5, N);

        cl_queue.enqueueNDRangeKernel(cl_intt, cl::NullRange, cl::NDRange(1), cl::NullRange);

        cl_queue.enqueueReadBuffer(cl_vals, CL_TRUE, 0, sizeof(clmp_uint64) * N, _vals);
    }

    void context::add_n_equal(clmpz_t *_lhs, clmpz_t *_rhs, clmp_uint64 _n_batch, clmpz_srcptr _mod) {
        cl::Buffer cl_lhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _n_batch);
        cl::Buffer cl_rhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _n_batch);
        cl::Buffer cl_mod(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t));

        cl_queue.enqueueWriteBuffer(cl_lhs, CL_TRUE, 0, sizeof(clmpz_t) * N * _n_batch, _lhs);
        cl_queue.enqueueWriteBuffer(cl_rhs, CL_TRUE, 0, sizeof(clmpz_t) * N * _n_batch, _rhs);
        cl_queue.enqueueWriteBuffer(cl_mod, CL_TRUE, 0, sizeof(clmpz_t), _mod);

        cl::Kernel cl_add_n_equal(cl_program, "add_n_equal");
        cl_add_n_equal.setArg(0, cl_lhs);
        cl_add_n_equal.setArg(1, cl_rhs);
        cl_add_n_equal.setArg(2, cl_mod);
        cl_add_n_equal.setArg(3, N);

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_add_n_equal, cl::NullRange, cl::NDRange(N, _n_batch), cl::NullRange);
        cl_queue.finish();
        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        std::cout << "add_n_equal kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_lhs, CL_TRUE, 0, sizeof(clmpz_t) * N * _n_batch, _lhs);
    }

    void context::right_shift_n_equal(clmpz_t *_vals, clmp_uint64 _n_batch, clmp_uint64 _bits) {
        cl::Buffer cl_vals(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _n_batch);
        cl::Buffer cl_tmp(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _n_batch);

        cl_queue.enqueueWriteBuffer(cl_vals, CL_TRUE, 0, sizeof(clmpz_t) * N * _n_batch, _vals);

        cl::Kernel cl_right_shift_n_equal(cl_program, "right_shift_n_equal");
        cl_right_shift_n_equal.setArg(0, cl_vals);
        cl_right_shift_n_equal.setArg(1, _bits);
        cl_right_shift_n_equal.setArg(2, N);
        cl_right_shift_n_equal.setArg(3, cl_tmp);

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_right_shift_n_equal, cl::NullRange, cl::NDRange(N, _n_batch), cl::NullRange);
        cl_queue.finish();
        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        std::cout << "right_shift_n_equal kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_vals, CL_TRUE, 0, sizeof(clmpz_t) * N * _n_batch, _vals);
    }

    void context::CRT(clmp_uint64 *_res, clmpz_t *_vals, const clmp_uint64 _n_batch, const long _np) {
        cl::Buffer cl_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * (_np << logN) * _n_batch);
        cl::Buffer cl_vals(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _n_batch);

        cl_queue.enqueueWriteBuffer(cl_vals, CL_TRUE, 0, sizeof(clmpz_t) * N * _n_batch, _vals);

        cl::Kernel cl_decomposition(cl_program, "crt_decomposition");
        cl_decomposition.setArg(0, cl_res);
        cl_decomposition.setArg(1, cl_vals);
        cl_decomposition.setArg(2, _np);
        cl_decomposition.setArg(3, logN);
        cl_decomposition.setArg(4, N);
        cl_decomposition.setArg(5, cl_p_vec);
        cl_decomposition.setArg(6, cl_p_r_vec);
        cl_decomposition.setArg(7, kbar2);

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_decomposition, cl::NullRange, cl::NDRange( N, _n_batch), cl::NullRange);
        cl_queue.finish();
        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        std::cout << "crt_decomposition kernel time : " << sec.count() << " seconds" << std::endl;

        cl::Kernel cl_crt_ntt(cl_program, "crt_ntt");
        cl_crt_ntt.setArg(0, cl_res);
        cl_crt_ntt.setArg(1, logN);
        cl_crt_ntt.setArg(2, N);
        cl_crt_ntt.setArg(3, cl_p_vec);
        cl_crt_ntt.setArg(4, cl_p_inv_vec);
        cl_crt_ntt.setArg(5, kbar2);
        cl_crt_ntt.setArg(6, cl_scaled_root_pows);

        start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_crt_ntt, cl::NullRange, cl::NDRange(_np), cl::NullRange);
        cl_queue.finish();
        sec = std::chrono::system_clock::now() - start;
        std::cout << "crt_ntt kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_res, CL_TRUE, 0, sizeof(clmp_uint64) * (_np << logN), _res);
    }

    void context::mult(clmpz_t *_res, const clmpz_t *_lhs, const clmpz_t *_rhs, const long _np, clmpz_srcptr _mod) {
        clmp_uint64 *r_lhs = new clmp_uint64[_np << logN];
        clmp_uint64 *r_rhs = new clmp_uint64[_np << logN];
        clmp_uint64 *r_res = new clmp_uint64[_np << logN];

        cl::Buffer cl_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N);
        cl::Buffer cl_lhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N);
        cl::Buffer cl_rhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N);
        cl::Buffer cl_mod(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t));

        cl::Buffer cl_tmp_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * (_np << logN));
        cl::Buffer cl_tmp_lhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * (_np << logN));
        cl::Buffer cl_tmp_rhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * (_np << logN));

        cl_queue.enqueueWriteBuffer(cl_lhs, CL_TRUE, 0, sizeof(clmpz_t) * N, _lhs);
        cl_queue.enqueueWriteBuffer(cl_rhs, CL_TRUE, 0, sizeof(clmpz_t) * N, _rhs);
        cl_queue.enqueueWriteBuffer(cl_mod, CL_TRUE, 0, sizeof(clmpz_t), _mod);

        cl::Kernel cl_decomposition(cl_program, "mult_decomposition");
        cl_decomposition.setArg(0, cl_tmp_lhs);
        cl_decomposition.setArg(1, cl_tmp_rhs);
        cl_decomposition.setArg(2, cl_lhs);
        cl_decomposition.setArg(3, cl_rhs);
        cl_decomposition.setArg(4, logN);
        cl_decomposition.setArg(5, cl_p_vec);
        cl_decomposition.setArg(6, cl_p_r_vec);
        cl_decomposition.setArg(7, kbar2);

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_decomposition, cl::NullRange, cl::NDRange(N, _np), cl::NullRange);
        cl_queue.finish();
        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        std::cout << "mult_decomposition kernel time : " << sec.count() << " seconds" << std::endl;

        cl::Kernel cl_ctx_mult(cl_program, "ctx_mult");
        cl_ctx_mult.setArg(0, cl_tmp_res);
        cl_ctx_mult.setArg(1, cl_tmp_lhs);
        cl_ctx_mult.setArg(2, cl_tmp_rhs);
        cl_ctx_mult.setArg(3, logN);
        cl_ctx_mult.setArg(4, N);
        cl_ctx_mult.setArg(5, cl_p_vec);
        cl_ctx_mult.setArg(6, cl_p_inv_vec);
        cl_ctx_mult.setArg(7, cl_p_r_vec);
        cl_ctx_mult.setArg(8, kbar2);
        cl_ctx_mult.setArg(9, cl_scaled_root_pows);
        cl_ctx_mult.setArg(10, cl_scaled_root_inv_pows);
        cl_ctx_mult.setArg(11, cl_scaled_N_inv);

        start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_ctx_mult, cl::NullRange, cl::NDRange(_np), cl::NullRange);
        cl_queue.finish();
        sec = std::chrono::system_clock::now() - start;
        std::cout << "ctx_mult kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_tmp_lhs, CL_TRUE, 0, sizeof(clmp_uint64) * (_np << logN), r_lhs);
        cl_queue.enqueueReadBuffer(cl_tmp_rhs, CL_TRUE, 0, sizeof(clmp_uint64) * (_np << logN), r_rhs);
        /*
        for (clmp_size_t p = 0; p < _np; p++) {
            for (clmp_uint64 i = 0; i < N; i++) {
                std::cout << r_lhs[p * N + i] << std::endl;
            }
        }
        std::cout << std::endl;

        for (clmp_size_t p = 0; p < _np; p++) {
            for (clmp_uint64 i = 0; i < N; i++) {
                std::cout << r_rhs[p * N + i] << std::endl;
            }
        }
        std::cout << std::endl << std::endl;
*/
        cl_queue.enqueueReadBuffer(cl_tmp_res, CL_TRUE, 0, sizeof(clmp_uint64) * (_np << logN), r_res);
/*
        for (clmp_size_t p = 0; p < _np; p++) {
            for (clmp_uint64 i = 0; i < N; i++) {
                std::cout << r_res[p * N + i] << std::endl;
            }
        }
        std::cout << std::endl;
*/
        cl::Kernel cl_reconstruct(cl_program, "reconstruct");
        cl_reconstruct.setArg(0, cl_res);
        cl_reconstruct.setArg(1, cl_tmp_res);
        cl_reconstruct.setArg(2, _np);
        cl_reconstruct.setArg(3, cl_mod);
        cl_reconstruct.setArg(4, nprimes);
        cl_reconstruct.setArg(5, cl_p_hat);
        cl_reconstruct.setArg(6, cl_p_hat_inv_mod_p);
        cl_reconstruct.setArg(7, cl_p_prod);
        cl_reconstruct.setArg(8, cl_p_prod_h);
        cl_reconstruct.setArg(9, cl_p_vec);
        cl_reconstruct.setArg(10, cl_p_r_vec);
        cl_reconstruct.setArg(11, kbar2);
        cl_reconstruct.setArg(12, N);
        cl_reconstruct.setArg(13, logN);

        start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_reconstruct, cl::NullRange, cl::NDRange(N), cl::NullRange);
        cl_queue.finish();
        sec = std::chrono::system_clock::now() - start;
        std::cout << "reconstruct kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_res, CL_TRUE, 0, sizeof(clmpz_t) * N, _res);
    }

    void context::mult_decrypt(clmpz_t *_res, const clmpz_t *_lhs, const clmpz_t *_rhs, const long _np, clmp_uint64 _n_batch, clmpz_srcptr _mod) {
        cl::Buffer cl_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _n_batch);
        cl::Buffer cl_lhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _n_batch);
        cl::Buffer cl_rhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N);
        cl::Buffer cl_mod(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t));

        cl::Buffer cl_tmp_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * (_np << logN) * _n_batch);
        cl::Buffer cl_tmp_lhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * (_np << logN) * _n_batch);
        cl::Buffer cl_tmp_rhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * (_np << logN));

        cl_queue.enqueueWriteBuffer(cl_lhs, CL_TRUE, 0, sizeof(clmpz_t) * N, _lhs);
        cl_queue.enqueueWriteBuffer(cl_rhs, CL_TRUE, 0, sizeof(clmpz_t) * N, _rhs);
        cl_queue.enqueueWriteBuffer(cl_mod, CL_TRUE, 0, sizeof(clmpz_t), _mod);

        cl::Kernel cl_decomposition(cl_program, "crt_decomposition");
        cl_decomposition.setArg(0, cl_tmp_lhs);
        cl_decomposition.setArg(1, cl_lhs);
        cl_decomposition.setArg(2, _np);
        cl_decomposition.setArg(3, logN);
        cl_decomposition.setArg(4, N);
        cl_decomposition.setArg(5, cl_p_vec);
        cl_decomposition.setArg(6, cl_p_r_vec);
        cl_decomposition.setArg(7, kbar2);

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_decomposition, cl::NullRange, cl::NDRange(N, _n_batch), cl::NullRange);
        cl_queue.finish();
        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        std::cout << "crt_decomposition kernel time : " << sec.count() << " seconds" << std::endl;

        cl_decomposition.setArg(0, cl_tmp_rhs);
        cl_decomposition.setArg(1, cl_rhs);
        cl_decomposition.setArg(2, _np);
        cl_decomposition.setArg(3, logN);
        cl_decomposition.setArg(4, N);
        cl_decomposition.setArg(5, cl_p_vec);
        cl_decomposition.setArg(6, cl_p_r_vec);
        cl_decomposition.setArg(7, kbar2);

        start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_decomposition, cl::NullRange, cl::NDRange(N, 1), cl::NullRange);
        cl_queue.finish();
        sec = std::chrono::system_clock::now() - start;
        std::cout << "crt_decomposition kernel time : " << sec.count() << " seconds" << std::endl;

        cl::Kernel cl_ctx_mult(cl_program, "ctx_mult_decrypt");
        cl_ctx_mult.setArg(0, cl_tmp_res);
        cl_ctx_mult.setArg(1, cl_tmp_lhs);
        cl_ctx_mult.setArg(2, cl_tmp_rhs);
        cl_ctx_mult.setArg(3, _np);
        cl_ctx_mult.setArg(4, logN);
        cl_ctx_mult.setArg(5, N);
        cl_ctx_mult.setArg(6, cl_p_vec);
        cl_ctx_mult.setArg(7, cl_p_inv_vec);
        cl_ctx_mult.setArg(8, cl_p_r_vec);
        cl_ctx_mult.setArg(9, kbar2);
        cl_ctx_mult.setArg(10, cl_scaled_root_pows);
        cl_ctx_mult.setArg(11, cl_scaled_root_inv_pows);
        cl_ctx_mult.setArg(12, cl_scaled_N_inv);

        start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_ctx_mult, cl::NullRange, cl::NDRange(_np, _n_batch), cl::NullRange);
        cl_queue.finish();
        sec = std::chrono::system_clock::now() - start;
        std::cout << "ctx_mult_decrypt kernel time : " << sec.count() << " seconds" << std::endl;

        cl::Kernel cl_reconstruct(cl_program, "reconstruct_batch");
        cl_reconstruct.setArg(0, cl_res);
        cl_reconstruct.setArg(1, cl_tmp_res);
        cl_reconstruct.setArg(2, _np);
        cl_reconstruct.setArg(3, cl_mod);
        cl_reconstruct.setArg(4, nprimes);
        cl_reconstruct.setArg(5, cl_p_hat);
        cl_reconstruct.setArg(6, cl_p_hat_inv_mod_p);
        cl_reconstruct.setArg(7, cl_p_prod);
        cl_reconstruct.setArg(8, cl_p_prod_h);
        cl_reconstruct.setArg(9, cl_p_vec);
        cl_reconstruct.setArg(10, cl_p_r_vec);
        cl_reconstruct.setArg(11, kbar2);
        cl_reconstruct.setArg(12, N);
        cl_reconstruct.setArg(13, logN);

        start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_reconstruct, cl::NullRange, cl::NDRange(N, _n_batch), cl::NullRange);
        cl_queue.finish();
        sec = std::chrono::system_clock::now() - start;
        std::cout << "reconstruct_batch kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_res, CL_TRUE, 0, sizeof(clmpz_t) * N * _n_batch, _res);
    }

    void context::mult_ntt(clmpz_t *_res, clmpz_t *_lhs, clmp_uint64 *_rhs, const clmp_uint64 _n_batch, const long _np, clmpz_srcptr _mod) {
        clmp_uint64 *r_lhs = new clmp_uint64[_np << logN];

        cl::Buffer cl_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _n_batch);
        cl::Buffer cl_lhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _n_batch);
        cl::Buffer cl_rhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * (_np << logN));
        cl::Buffer cl_mod(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t));

        cl::Buffer cl_tmp_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * (_np << logN) * _n_batch);
        cl::Buffer cl_tmp_lhs(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * (_np << logN) * _n_batch);

        cl_queue.enqueueWriteBuffer(cl_lhs, CL_TRUE, 0, sizeof(clmpz_t) * N * _n_batch, _lhs);
        cl_queue.enqueueWriteBuffer(cl_rhs, CL_TRUE, 0, sizeof(clmp_uint64) * (_np << logN) * _n_batch, _rhs);
        cl_queue.enqueueWriteBuffer(cl_mod, CL_TRUE, 0, sizeof(clmpz_t), _mod);

        cl::Kernel cl_decomposition(cl_program, "mult_ntt_decomposition");
        cl_decomposition.setArg(0, cl_tmp_res);
        cl_decomposition.setArg(1, cl_lhs);
        cl_decomposition.setArg(2, _np);
        cl_decomposition.setArg(3, logN);
        cl_decomposition.setArg(4, N);
        cl_decomposition.setArg(5, cl_p_vec);
        cl_decomposition.setArg(6, cl_p_r_vec);
        cl_decomposition.setArg(7, kbar2);

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_decomposition, cl::NullRange, cl::NDRange(N, _n_batch), cl::NullRange);
        cl_queue.finish();
        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        std::cout << "crt_decomposition kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_tmp_res, CL_TRUE, 0, sizeof(clmp_uint64) * (_np << logN), r_lhs);

        for (clmp_uint64 i = 0; i < (_np << logN); i++) {
            std::cout << r_lhs[i] << std::endl;
        }
        std::cout << std::endl;

        cl::Kernel cl_ctx_mult(cl_program, "ctx_mult_partial_ntt");
        cl_ctx_mult.setArg(0, cl_tmp_res);
        cl_ctx_mult.setArg(1, cl_tmp_lhs);
        cl_ctx_mult.setArg(2, cl_rhs);
        cl_ctx_mult.setArg(3, _np);
        cl_ctx_mult.setArg(4, logN);
        cl_ctx_mult.setArg(5, N);
        cl_ctx_mult.setArg(6, cl_p_vec);
        cl_ctx_mult.setArg(7, cl_p_inv_vec);
        cl_ctx_mult.setArg(8, cl_p_r_vec);
        cl_ctx_mult.setArg(9, kbar2);
        cl_ctx_mult.setArg(10, cl_scaled_root_pows);
        cl_ctx_mult.setArg(11, cl_scaled_root_inv_pows);
        cl_ctx_mult.setArg(12, cl_scaled_N_inv);

        start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_ctx_mult, cl::NullRange, cl::NDRange(_np, _n_batch), cl::NullRange);
        cl_queue.finish();
        sec = std::chrono::system_clock::now() - start;
        std::cout << "ctx_mult_partial_ntt kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_tmp_res, CL_TRUE, 0, sizeof(clmp_uint64) * (_np << logN), r_lhs);

        for (clmp_uint64 i = 0; i < (_np << logN); i++) {
            std::cout << r_lhs[i] << std::endl;
        }
        std::cout << std::endl;

        cl::Kernel cl_reconstruct(cl_program, "reconstruct_batch");
        cl_reconstruct.setArg(0, cl_res);
        cl_reconstruct.setArg(1, cl_tmp_res);
        cl_reconstruct.setArg(2, _np);
        cl_reconstruct.setArg(3, cl_mod);
        cl_reconstruct.setArg(4, nprimes);
        cl_reconstruct.setArg(5, cl_p_hat);
        cl_reconstruct.setArg(6, cl_p_hat_inv_mod_p);
        cl_reconstruct.setArg(7, cl_p_prod);
        cl_reconstruct.setArg(8, cl_p_prod_h);
        cl_reconstruct.setArg(9, cl_p_vec);
        cl_reconstruct.setArg(10, cl_p_r_vec);
        cl_reconstruct.setArg(11, kbar2);
        cl_reconstruct.setArg(12, N);
        cl_reconstruct.setArg(13, logN);

        start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_reconstruct, cl::NullRange, cl::NDRange(N, _n_batch), cl::NullRange);
        cl_queue.finish();
        sec = std::chrono::system_clock::now() - start;
        std::cout << "reconstruct kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_res, CL_TRUE, 0, sizeof(clmpz_t) * N * _n_batch, _res);
    }

    void context::test(clmpz_t * x, clmpz_t * a, clmpz_t * b, long np, clmpz_srcptr mod) {
        clmp_uint64 *rx = new clmp_uint64[np << logN];
        clmp_uint64 *ax = new clmp_uint64[np << logN];
        clmp_uint64 *bx = new clmp_uint64[np << logN];

        for (long i = 0; i < np; i++) {
            clmp_uint64 *rxi = rx + (i << logN);
            clmp_uint64 *rai = ax + (i << logN);
            clmp_uint64 *rbi = bx + (i << logN);
            clmp_uint64 pi = p_vec[i];

            for (clmp_uint64 n = 0; n < N; n++) {
                clmpz_t tmp;
                clmpz_mod_ul(tmp, a[n], pi);
                rai[n] = clmpz_to_ul(tmp);
                clmpz_mod_ul(tmp, b[n], pi);
                rbi[n] = clmpz_to_ul(tmp);
            }

            NTT_cpu(rai, i);
            NTT_cpu(rbi, i);

            for (clmp_uint64 n = 0; n < N; n++) {
                rxi[n] = mul_mod(rai[n], rbi[n], pi);
            }
            iNTT_cpu(rxi, i);

        }

        reconstruct_cpu(x, rx, np, mod);
        /*
        cl::Buffer cl_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N);
        cl::Buffer cl_vals(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmp_uint64) * (np << logN));
        cl::Buffer cl_mod(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t));

        cl_queue.enqueueWriteBuffer(cl_vals, CL_TRUE, 0, sizeof(clmp_uint64) * (np << logN), rx);
        cl_queue.enqueueWriteBuffer(cl_mod, CL_TRUE, 0, sizeof(clmpz_t), mod);


        cl::Kernel cl_reconstruct(cl_program, "reconstruct");
        cl_reconstruct.setArg(0, cl_res);
        cl_reconstruct.setArg(1, cl_vals);
        cl_reconstruct.setArg(2, np);
        cl_reconstruct.setArg(3, cl_mod);
        cl_reconstruct.setArg(4, nprimes);
        cl_reconstruct.setArg(5, cl_p_hat);
        cl_reconstruct.setArg(6, cl_p_hat_inv_mod_p);
        cl_reconstruct.setArg(7, cl_p_prod);
        cl_reconstruct.setArg(8, cl_p_prod_h);
        cl_reconstruct.setArg(9, cl_p_vec);
        cl_reconstruct.setArg(10, cl_p_r_vec);
        cl_reconstruct.setArg(11, kbar2);
        cl_reconstruct.setArg(12, N);
        cl_reconstruct.setArg(13, logN);

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_reconstruct, cl::NullRange, cl::NDRange(1), cl::NullRange);
        cl_queue.finish();
        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        std::cout << "kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_res, CL_TRUE, 0, sizeof(clmpz_t) * N, x);
*/
        delete [] rx;
        delete [] ax;
        delete [] bx;

    }

    void context::reconstruct_cpu(clmpz_t *_res, clmp_uint64* _vals, const long _np, clmpz_srcptr _q) {
        clmpz_t *p_hat_np = p_hat + (_np - 1) * nprimes;
        clmp_uint64 *p_hat_inv_mod_p_np = p_hat_inv_mod_p + (_np - 1) * nprimes;
        clmpz_t *p_prod_np = p_prod + _np - 1;
        clmpz_t *p_prod_h_np = p_prod_h + _np - 1;

        for (clmp_uint64 n = 0; n < N; n++) {
            clmpz_set_ui(_res[n], 0);

            for (long i = 0; i < _np; i++) {
                clmpz_t tmp, tmp2;
                clmp_uint64 s = mul_mod(_vals[n + (i << logN)], *(p_hat_inv_mod_p_np + i), p_vec[i]);
                clmpz_mul_ul(tmp, *(p_hat_np + i), s);
                clmpz_set(tmp2, _res[n]);
                clmpz_add(_res[n], tmp2, tmp);
            }

            clmpz_mod(_res[n], _res[n], *p_prod_np);

            if (clmpz_cmp(_res[n], *p_prod_h_np) == 1) {
                clmpz_sub(_res[n], _res[n], *p_prod_np);
            }

            clmpz_mod(_res[n], _res[n], _q);
        }
    }

    void context::NTT_cpu(clmp_uint64 *_vals, long index) {
        long t = N;
        long logt1 = logN + 1;
        clmp_uint64 p = p_vec[index];
        clmp_uint64 pInv = p_inv_vec[index];

        for (clmp_uint64 m = 1; m < N; m <<= 1) {
            t >>= 1;
            logt1 -= 1;
            for (clmp_uint64 i = 0; i < m; i++) {
                long j1 = i << logt1;
                long j2 = j1 + t - 1;
                clmp_uint64 W = scaled_root_pows[index * N + m + i];
                for (long j = j1; j <= j2; j++) {
                    butt(_vals[j], _vals[j+t], W, p, pInv);
                }
            }
        }
    }

    void context::iNTT_cpu(clmp_uint64 *_vals, long index) {
        clmp_uint64 p = p_vec[index];
        clmp_uint64 pInv = p_inv_vec[index];
        long t = 1;
        for (clmp_uint64 m = N; m > 1; m >>= 1) {
            long j1 = 0;
            long h = m >> 1;
            for (long i = 0; i < h; i++) {
                long j2 = j1 + t - 1;
                clmp_uint64 W = scaled_root_inv_pows[index * N + h + i];
                for (long j = j1; j <= j2; j++) {
                    ibutt(_vals[j], _vals[j+t], W, p, pInv);
                }
                j1 += (t << 1);
            }
            t <<= 1;
        }

        clmp_uint64 NScale = scaled_N_inv[index];
        for (clmp_uint64 i = 0; i < N; i++) {
            idivN(_vals[i], NScale, p, pInv);
        }
    }

    void context::butt(clmp_uint64& a, clmp_uint64& b, clmp_uint64 W, clmp_uint64 p, clmp_uint64 pInv) {
        unsigned __int128 U = static_cast<unsigned __int128>(b) * W;
        clmp_uint64 U0 = static_cast<clmp_uint64>(U);
        clmp_uint64 U1 = U >> 64;
        clmp_uint64 Q = U0 * pInv;
        unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
        clmp_uint64 H = Hx >> 64;
        clmp_uint64 V = U1 < H ? U1 + p - H : U1 - H;
        b = a < V ? a + p - V : a - V;
        a += V;
        if (a > p) a -= p;
    }

    void context::ibutt(clmp_uint64& a, clmp_uint64& b, clmp_uint64 W, clmp_uint64 p, clmp_uint64 pInv) {
        clmp_uint64 T = a < b ? a + p - b : a - b;
        a += b;
        if (a > p) a -= p;
        unsigned __int128 UU = static_cast<unsigned __int128>(T) * W;
        clmp_uint64 U0 = static_cast<clmp_uint64>(UU);
        clmp_uint64 U1 = UU >> 64;
        clmp_uint64 Q = U0 * pInv;
        unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
        clmp_uint64 H = Hx >> 64;
        b = (U1 < H) ? U1 + p - H : U1 - H;
    }

    void context::idivN(clmp_uint64& a, clmp_uint64 NScale, clmp_uint64 p, clmp_uint64 pInv) {
        unsigned __int128 U = static_cast<unsigned __int128>(a) * NScale;
        clmp_uint64 U0 = static_cast<clmp_uint64>(U);
        clmp_uint64 U1 = U >> 64;
        clmp_uint64 Q = U0 * pInv;
        unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
        clmp_uint64 H = Hx >> 64;
        a = (U1 < H) ? U1 + p - H : U1 - H;
    }

    void context::sub_from_gauss_n_equal(clmpz_t *_res, clmpz_srcptr _mod) {
        clmp_uint64 seed = clmp_rand_bnd_ui(0xFFFFFFF);

        cl::Buffer cl_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N);
        cl::Buffer cl_mod(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t));

        cl_queue.enqueueWriteBuffer(cl_res, CL_TRUE, 0, sizeof(clmpz_t) * N, _res);
        cl_queue.enqueueWriteBuffer(cl_mod, CL_TRUE, 0, sizeof(clmpz_t), _mod);

        cl::Kernel cl_sub_from_gauss(cl_program, "sub_from_gauss_n_equal");
        cl_sub_from_gauss.setArg(0, cl_res);
        cl_sub_from_gauss.setArg(1, cl_mod);
        cl_sub_from_gauss.setArg(2, seed);
        cl_sub_from_gauss.setArg(3, sigma);

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_sub_from_gauss, cl::NullRange, cl::NDRange(N), cl::NullRange);
        cl_queue.finish();
        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        std::cout << "sub_from_gauss_n_equal kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_res, CL_TRUE, 0, sizeof(clmpz_t) * N, _res);
    }

    void context::sub_from_gauss_n_equal_cpu(clmpz_t *_res, clmpz_srcptr _mod) {
        static clmp_int32 bignum = 0xFFFFFFF;

        for (clmp_uint64 i = 0; i < N; i += 2) {
            double r1 = (1 + clmp_rand_bnd_ui(bignum)) / (__CLMP_CAST(double, bignum) + 1);
            double r2 = (1 + clmp_rand_bnd_ui(bignum)) / (__CLMP_CAST(double, bignum) + 1);
            double theta = 2 * M_PI * r1;
            double rr = sqrt(-2.0 * log(r2)) * sigma;

            clmpz_t tmp;
            clmpz_neg(tmp, _res[i]);
            clmpz_add_ul(_res[i], tmp, __CLMP_CAST(clmp_int64, floor(rr * cos(theta) + 0.5)));
            clmpz_mod(_res[i], _res[i], _mod);

            clmpz_neg(tmp, _res[i + 1]);
            clmpz_add_ul(_res[i + 1], tmp, __CLMP_CAST(clmp_int64, floor(rr * sin(theta) + 0.5)));
            clmpz_mod(_res[i + 1], _res[i + 1], _mod);
        }
    }

    void context::add_gauss_n_equal(clmpz_t *_res, clmp_uint64 _n_batch, clmpz_srcptr _mod) {
        clmp_uint64 seed = clmp_rand_bnd_ui(0xFFFFFFF);

        cl::Buffer cl_res(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t) * N * _n_batch);
        cl::Buffer cl_mod(*cl_ctx, CL_MEM_READ_WRITE, sizeof(clmpz_t));

        cl_queue.enqueueWriteBuffer(cl_res, CL_TRUE, 0, sizeof(clmpz_t) * N * _n_batch, _res);
        cl_queue.enqueueWriteBuffer(cl_mod, CL_TRUE, 0, sizeof(clmpz_t), _mod);

        cl::Kernel cl_add_gauss(cl_program, "add_gauss_n_equal");
        cl_add_gauss.setArg(0, cl_res);
        cl_add_gauss.setArg(1, cl_mod);
        cl_add_gauss.setArg(2, seed);
        cl_add_gauss.setArg(3, sigma);
        cl_add_gauss.setArg(4, N);

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
        cl_queue.enqueueNDRangeKernel(cl_add_gauss, cl::NullRange, cl::NDRange(N, _n_batch), cl::NullRange);
        cl_queue.finish();
        std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
        std::cout << "add_gauss_n_equal kernel time : " << sec.count() << " seconds" << std::endl;

        cl_queue.enqueueReadBuffer(cl_res, CL_TRUE, 0, sizeof(clmpz_t) * N * _n_batch, _res);
    }

    void context::add_gauss_n_equal_cpu(clmpz_t *_res, clmpz_srcptr _mod) {
        static clmp_int32 bignum = 0xFFFFFFF;

        for (clmp_uint64 i = 0; i < N; i += 2) {
            double r1 = (1 + clmp_rand_bnd_ui(bignum)) / (__CLMP_CAST(double, bignum) + 1);
            double r2 = (1 + clmp_rand_bnd_ui(bignum)) / (__CLMP_CAST(double, bignum) + 1);
            double theta = 2 * M_PI * r1;
            double rr = sqrt(-2.0 * log(r2)) * sigma;

            clmpz_add_ul(_res[i], _res[i], __CLMP_CAST(clmp_int64, floor(rr * cos(theta) + 0.5)));
            clmpz_mod(_res[i], _res[i], _mod);

            clmpz_add_ul(_res[i + 1], _res[i + 1], __CLMP_CAST(clmp_int64, floor(rr * sin(theta) + 0.5)));
            clmpz_mod(_res[i + 1], _res[i + 1], _mod);
        }
    }

    void context::sampleHWT(clmpz_t *_res) {
        clmp_uint64 idx = 0;
        clmpz_t tmp;
        clmpz_rand_bits(tmp, h);
        while (idx < h) {
            long i = clmp_rand_bits_ui(logN);

            if (clmpz_is_zero(_res[i])) {
                if (clmpz_bit(tmp, idx) == 0) {;
                    clmpz_set_i(_res[i], 1);
                }
                else {
                    clmpz_set_i(_res[i], -1);
                }
                idx++;
            }
        }
    }

    void context::sample_ZO(clmpz_t *_res) {
        clmpz_t tmp;

        clmpz_rand_bits(tmp, M);

        for (clmp_uint64 i = 0; i < N; i++) {
            if (clmpz_bit(tmp, 2 * i) == 0) {
                clmpz_set_i(_res[i], 0);
            }
            else {
                if (clmpz_bit(tmp, 2 * i + 1) == 0) {
                    clmpz_set_i(_res[i], 1);
                }
                else {
                    clmpz_set_i(_res[i], -1);
                }
            }
        }
    }

    void context::sample_uniform(clmpz_t *_res, long _bits) {
        for (clmp_uint64 i = 0; i < N; i++) {
            clmpz_rand_bits(_res[i], _bits);
        }
    }
}
