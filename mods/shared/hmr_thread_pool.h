#ifndef HMR_THREAD_POOL_H
#define HMR_THREAD_POOL_H

#include <thread>
#include <atomic>
#include <functional>
#include <condition_variable>

namespace hmr
{
    typedef std::int_fast32_t   hmr_i32;
    typedef std::uint_fast32_t  hmr_ui32;
    typedef std::uint_fast64_t  hmr_ui64;

    template <typename T>
    class thread_pool_queue
    {
    public:
        thread_pool_queue(const hmr_i32 size_limit) :
            m_queue(new T[size_limit]),
            m_queueSize(size_limit),
            m_head(0),
            m_tail(0)
        {
        }

        T& front()
        {
            return m_queue[m_head];
        }

        const T& front() const
        {
            return m_queue[m_head];
        }

        bool empty() const
        {
            return m_head == m_tail;
        }

        bool full() const
        {
            return (m_tail + 1 == m_head) || ((m_tail == 0) && (m_head == m_queueSize - 1));
        }

        void pop()
        {
            //Only pop when 
            if (empty())
            {
                throw std::runtime_error("Thread pool queue is empty to pop a task.");
            }
            //Move the head to next position.
            m_head = next_pos(m_head);
        }

        void push(const T& element)
        {
            if (full())
            {
                throw std::runtime_error("Thread pool queue is already full to push a task.");
            }
            //Copy the data to the tail.
            m_queue[m_tail] = element;
            //Move the tail to next position.
            m_tail = next_pos(m_tail);
        }

        ~thread_pool_queue()
        {
            delete[] m_queue;
        }

    private:
        inline size_t next_pos(size_t current)
        {
            ++current;
            return current >= m_queueSize ? 0 : current;
        }
        T* m_queue;
        size_t m_head, m_tail, m_queueSize;
    };

    template <typename T>
    class thread_pool
    {
    public:
        thread_pool(void (*task)(const T&), const hmr_i32 queue_limit, const hmr_ui32 threads = std::thread::hardware_concurrency()) :
            m_taskRequests(queue_limit),
            m_task(task),
            m_threadCount(threads)
        {
            //Create threads.
            m_threads = new std::thread[m_threadCount];
            for (hmr_ui32 i = 0; i < m_threadCount; ++i)
            {
                m_threads[i] = std::thread(&thread_pool::worker, this);
            }
        }

        ~thread_pool()
        {
            //Wait for all the tasks complete.
            wait_for_tasks();
            //Shutdown the pool.
            m_running = false;
            //Destroy all the tasks.
            for (hmr_ui32 i = 0; i < m_threadCount; ++i)
            {
                m_threads[i].join();
            }
            delete[] m_threads;
        }

        void push_task(const T& task)
        {
            //Wait for the queue is valid.
            std::unique_lock<std::mutex> queue_lock(m_queueMutex);
            m_queuePushCv.wait(queue_lock, [this]{ return !m_taskRequests.full(); });
            //Push the data back.
            ++m_tasksTotal;
            m_taskRequests.push(task);
        }

        void wait_for_tasks()
        {
            while (true)
            {
                //When the tasks are completed, then done.
                if (m_tasksTotal == 0)
                {
                    break;
                }
                //Hold and wait for a while.
                handover_core();
            }
        }

    private:
        bool assign_task(T& task)
        {
            const std::unique_lock<std::mutex> queue_lock(m_queueMutex);
            //Try to pop the task.
            if (m_taskRequests.empty())
            {
                return false;
            }
            //Assign the task data.
            task = std::move(m_taskRequests.front());
            m_taskRequests.pop();
            m_queuePushCv.notify_one();
            return true;
        }

        void worker()
        {
            //While the thread pool is running, try to fetch the task parameter from the queue.
            while (m_running)
            {
                T task_parameter;
                if (assign_task(task_parameter))
                {
                    m_task(task_parameter);
                    --m_tasksTotal;
                    m_queueEmpty.notify_one();
                }
                else
                {
                    handover_core();
                }
            }
        }

        inline void handover_core()
        {
            //std::this_thread::yield();
            std::this_thread::sleep_for(std::chrono::microseconds(1000));
        }

        mutable std::mutex m_queueMutex = {};
        std::condition_variable m_queuePushCv, m_queueEmpty;

        thread_pool_queue<T> m_taskRequests;
        std::atomic<bool> m_running{ true };
        std::atomic<hmr_ui32> m_tasksTotal{ 0 };
        void (*m_task)(const T&);
        std::thread* m_threads;
        hmr_ui32 m_threadCount;
    };
}

#endif // HMR_THREAD_POOL_H
