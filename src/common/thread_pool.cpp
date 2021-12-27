#include "common/thread_pool.h"

namespace COMMON {
    thread_pool::thread_pool(size_t _num_thread):__num_thread(_num_thread), __stop_all(false) {
        __workers.reserve(__num_thread);

        for (size_t i = 0; i < __num_thread; i++) {
            __workers.emplace_back([this]() { this->run_thread(); });
        }
    }

    thread_pool::~thread_pool() {
        __stop_all = true;
        __cv_job_q.notify_all();

        for (auto &e: __workers) {
            e.join();
        }
    }

    template <class F, class... Args>
    std::future<typename std::result_of<F(Args...)>::type> enqueue_job(F&& f, Arg&&... args) {
        if (__stop_all) {
            throw std::runtime_error("threads are suspended.");
        }

        using return_type = typename std::result_of<F(Args...)>::type;
        auto job = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));
        std::future<return_type> job_result_future = job->get_future();
        {
            std::lock_guard<std::mutex> lock(__m_job_q_);
            jobs_.push([job]() { (*job)(); });
        }
        __cv_job_q_.notify_one();

        return job_result_future;
    }

    void thread_pool::run_thread() {
        while(true) {
            std::unique_lock<std::mutex> lock(__m_job_q);
            __cv_job_q.wait(lock, [this]() { return !this->__jobs.empty() || __stop_all; });

            if (__stop_all && this->__jobs.empty()) {
                return;
            }

            std::function<void()> job = std::move(__jobs.front());
            __jobs.pop();
            lock.unlock();

            job();
        }
    }
}
