#ifndef HMR_THREAD_POOL_H
#define HMR_THREAD_POOL_H

#include <atomic>
#include <list>
#include <chrono>
#include <cstdint>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>

template <typename F, typename T>
class thread_pool
{
    typedef std::uint_fast32_t ui32;
    typedef std::uint_fast64_t ui64;

public:
    thread_pool(F func, const ui32 &_thread_count = std::thread::hardware_concurrency()):
        func(func),
        thread_count(_thread_count ? _thread_count : std::thread::hardware_concurrency()), threads(new std::thread[_thread_count ? _thread_count : std::thread::hardware_concurrency()])
    {
        create_threads();
    }

    ~thread_pool()
    {
        wait_for_tasks();
        running = false;
        destroy_threads();
    }

    void push_task(const T &task)
    {
        tasks_total++;
        {
            const std::unique_lock<std::mutex> lock(queue_mutex);
            tasks.push(task);
        }
    }

    void wait_for_tasks()
    {
        while (true)
        {
            if (!paused)
            {
                if (tasks_total == 0)
                    break;
            }
            else
            {
                if (get_tasks_running() == 0)
                    break;
            }
            sleep_or_yield();
        }
    }

    std::atomic<bool> paused {false};
    ui32 sleep_duration = 1000;

private:
    ui32 get_tasks_running() const
    {
        return tasks_total - (ui32)get_tasks_queued();
    }

    ui64 get_tasks_queued() const
    {
        const std::unique_lock<std::mutex> lock(queue_mutex);
        return tasks.size();
    }

    virtual void create_threads()
    {
        for (ui32 i = 0; i < thread_count; i++)
        {
            threads[i] = std::thread(&thread_pool::worker, this);
        }
    }

    void destroy_threads()
    {
        for (ui32 i = 0; i < thread_count; i++)
        {
            threads[i].join();
        }
    }

    bool pop_task(T &task)
    {
        const std::unique_lock<std::mutex> lock(queue_mutex);
        if (tasks.empty())
            return false;
        else
        {
            task = std::move(tasks.front());
            tasks.pop();
            return true;
        }
    }

    void sleep_or_yield()
    {
        if (sleep_duration)
            std::this_thread::sleep_for(std::chrono::microseconds(sleep_duration));
        else
            std::this_thread::yield();
    }

    void worker()
    {
        while (running)
        {
            T param;
            if (!paused && pop_task(param))
            {
                func(param);
                tasks_total--;
            }
            else
            {
                sleep_or_yield();
            }
        }
    }

    mutable std::mutex queue_mutex = {};
    F *func;
    std::atomic<bool> running {true};
    std::queue<T> tasks;
    ui32 thread_count;
    std::unique_ptr<std::thread[]> threads;
    std::atomic<ui32> tasks_total {0};
};

#endif // HMR_THREAD_POOL_H
