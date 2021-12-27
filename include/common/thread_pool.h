#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

namespace COMMON {
    class thread_pool {
    public:
        thread_pool(size_t _num_thread);
        ~thread_pool();

        template <class F, class... Args>
        std::future<typename std::result_of<F(Args...)>::type> enqueue_job(F&& f, Args&&... args);

    private:
        size_t __num_thread;

        std::vector<std::thread> __workers;

        std::queue<std::function<void()>> __jobs;

        std::condition_variable __cv_job_q;
        std::mutex __m_job_q;

        bool __stop_all;

        void run_thread();
    };
}


#endif
